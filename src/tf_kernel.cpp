#include <Rcpp.h>
#include <omp.h>
/*
 * using an R data.frame object, compute tfIdf/cosine kernel
 * Source code for calculating term frequencecy kernel using the given TF <-> gene data
 */

/*
 * store stuff
 */
class Tfdatamat{
private:
	std::vector<std::string> genes,tfs;
	SEXP num_mat;
public:
	Tfdatamat(std::vector<std::string> genes,std::vector<std::string> tfs,SEXP num_mat);
//	some getters
	std::vector<std::string> get_genes();
	std::vector<std::string> get_tfs();
	SEXP get_mat();
};

Tfdatamat::Tfdatamat(std::vector<std::string> genes,std::vector<std::string> tfs,SEXP num_mat){
	this->genes = genes;
	this->tfs = tfs;
	this->num_mat = num_mat;
}

std::vector<std::string> Tfdatamat::get_genes(){
	return genes;
}

std::vector<std::string> Tfdatamat::get_tfs(){
	return tfs;
}

SEXP Tfdatamat::get_mat(){
	return num_mat;
}

/*
 * You need these variables set to enable openmp
 * Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
 * Sys.setenv("PKG_LIBS"="-fopenmp")
 */

/*
 * remove duplicates
 * shamelessly copied from here:
 * http://stackoverflow.com/questions/1041620/whats-the-most-efficient-way-to-erase-duplicates-and-sort-a-vector
 */
std::vector<std::string> dedup(std::vector<std::string> in_v){
	std::set<std::string> nodup;
	std::vector<std::string> nodup_vec;
	unsigned size = in_v.size();
	for(unsigned i = 0; i < size; ++i){
		nodup.insert(in_v[i]);
	}
	nodup_vec.assign(nodup.begin(), nodup.end());
	return nodup_vec;
}

/*
 * given two vectors, get the count of common terms in them
 * modified from here:
 * @param v1: string vector
 * @param v2: another string vector
 * @param intersect: boolean, if true, find length of  common terms, if false return the length of deduplicated union
 * http://stackoverflow.com/questions/2404094/how-do-i-get-characters-common-to-two-vectors-in-c
 */
int getcount(std::vector<std::string> v1, std::vector<std::string> v2,bool intersect=true){
	int count=0;
	if(intersect){
		for(std::vector<std::string>::iterator m=v1.begin();m!=v1.end();++m){
				if(std::find(v2.begin(),v2.end(),*m)!=v2.end()){
					count++;
				}
			}
	}else{
		v1.insert(v1.end(),v2.begin(),v2.end());
		v1 = dedup(v1);
		count = v1.size();
	}

	return count;
}


/*
 * given a string vector, count number of times a string appears
 * another shameless  copy:
 * http://www.cplusplus.com/reference/algorithm/count/
 */
std::map<std::string,int> getcount(std::vector<std::string> vs){
	std::vector<std::string> vs_dedup = dedup(vs);
	std::map<std::string,int> v_map;
	for(int i=0;i<vs_dedup.size();i++){
		int count = std::count(vs.begin(),vs.end(),vs_dedup[i]);
		v_map[vs_dedup[i]]=count;
	}
	return v_map;
}

/*
 * given a map as input, return keys as arrays
 */
std::vector<std::string> getkeys(std::map<std::string,std::vector<std::string> > rmap){
	std::vector<std::string> keymap;
	std::map<std::string,std::vector<std::string> >::iterator it;
	for(it=rmap.begin();it!=rmap.end();++it){
		keymap.push_back(it->first);
	}
	return keymap;
}

/*
 * given a map of vectors, calculate term frequency
 * here map key is gene id and value string vector contains tfs per gene
 */
std::map<std::string,std::map<std::string,int> > get_tf(std::map<std::string,std::vector<std::string> > gtf_map){
//	return me
	std::map<std::string,std::map<std::string,int> > gtf_cmap;
	std::map<std::string,std::vector<std::string> >::iterator it;
	for( it=gtf_map.begin(); it!=gtf_map.end();++it){
		std::string gene = it->first;
		std::vector<std::string> tfs = it->second;
		std::map<std::string,int> tf_cmap = getcount(tfs);
		gtf_cmap[gene]=tf_cmap;
	}
	return gtf_cmap;
}

/*
 * given two vectors calculate inverse document frequency
 * first vector stands for document or gene and
 * the second vector stands for word or tf
 */
std::map<std::string,double> get_idf(std::vector<std::string> g,std::vector<std::string>tf){
	std::map<std::string,std::vector<std::string> >tfg_map;
//	map transcription factor to gene
	for(int i=0;i<tf.size();i++){
		tfg_map[tf[i]].push_back(g[i]);
	}
//	total documents(genes)
	int total_genes = dedup(g).size();
//	sort and count, to get the number of genes per tf AND return!!!!
	std::map<std::string,double> tfg_cmap;
	std::map<std::string,std::vector<std::string> >::iterator tgifit;
	for(tgifit=tfg_map.begin();tgifit!=tfg_map.end();++tgifit){
		std::string tf = tgifit->first;
		int gcount = dedup(tgifit->second).size();
		double idf = log(double(total_genes)/double(gcount));
		tfg_cmap[tf]= idf;
	}
	return tfg_cmap;
}

/*
 * use calculated term frequency and inverse document frequency to calculate tfidf
 * @param tfmap: term frequency map per document (gene), per word (transcription factor)
 * @param idfmap: inverse document frequency per word (transcription factor)
 * @param notf: boolean, whether to return matrix with inverse document frequencies only
 * @param binary: boolean, whether to return binary matrix ( value 1 if tf binds to gene, 0 if not)
 */

Tfdatamat calc_tfidf(std::map<std::string,std::map<std::string,int> > tfmap, std::map<std::string,double> idfmap,bool notf = false, bool binary = false){
//	return a matrix
//	define number of rows and cols first, which should be:
//	rows: the number of genes and (so tfmap, term frequency map size)
//	cols: the number of transcription factors (so idfmap, inverse document frequency map size)
	Rcpp::NumericMatrix tfidf_mat (tfmap.size(), idfmap.size());
//	iterator for tfmap
	std::map<std::string,std::map<std::string,int> >::iterator tfit;
//	iterator for idfmap
	std::map<std::string,double>::iterator idfit;
//	the row counter
	int row=0;
//	define vectors for gene and tf ids
	std::vector<std::string> genes,tfs;
//	loopy
	for(tfit=tfmap.begin();tfit!=tfmap.end();++tfit){
		std::string gene = tfit->first;
		genes.push_back(gene);
		std::map<std::string,int> tfcount_map = tfit->second;
//		the column counter
		int col=0;
		for(idfit=idfmap.begin();idfit!=idfmap.end();++idfit){
			std::string tf = idfit->first;
			if(row==0){
				tfs.push_back(tf); // add the tfs only for the first iteration
			}
			double idf = idfit->second;
			double tfidf = idf;
			std::map<std::string,int>::iterator tf_find = tfcount_map.find(tf);
			if(tf_find!=tfcount_map.end()){
				int tf = tf_find->second;
				if(binary){
					tfidf = 1.00;
				}else if(notf){
					tfidf = double(1) * idf;
				}else{
					tfidf = log(double(tf)+1) * idf;
				}
			}else{
				tfidf=0;
			}
			tfidf_mat(row,col)=tfidf;
			col++;
		}
		row++;
	}
	Tfdatamat tdd(genes,tfs,tfidf_mat);
	return tdd;
}


// [[Rcpp::export]]
SEXP get_tfidf(SEXP df,bool notf){
	Rcpp::DataFrame tf_df(df);
	Rcpp::StringVector g_(tf_df[0]);
	Rcpp::StringVector tf_(tf_df[1]);
	std::vector<std::string> g = Rcpp::as<std::vector<std::string> >(g_);
	std::vector<std::string> tf = Rcpp::as<std::vector<std::string> >(tf_);
//	gene to tf map
	std::map<std::string,std::vector<std::string> >gtf_map;
	for(int i=0;i<g.size();i++){
		gtf_map[g[i]].push_back(tf[i]);
	}
// now count
	std::map<std::string,std::map<std::string,int> >cmap = get_tf(gtf_map);
	std::map<std::string,double> tfg_idfmap = get_idf(g,tf);
//	return some stuff
	Tfdatamat calc_mat = calc_tfidf(cmap,tfg_idfmap,notf);
	Rcpp::List out_data = Rcpp::List::create(
			Rcpp::Named("mat")=calc_mat.get_mat(),
			Rcpp::Named("genes")=calc_mat.get_genes(),
			Rcpp::Named("tfs")=calc_mat.get_tfs()
			);
	return Rcpp::wrap(out_data);
}

/*
 * get binary matrix
 */

// [[Rcpp::export]]
SEXP get_tf_gene_binary(SEXP df){
	Rcpp::DataFrame tf_df(df);
	Rcpp::StringVector g_(tf_df[0]);
	Rcpp::StringVector tf_(tf_df[1]);
	std::vector<std::string> g = Rcpp::as<std::vector<std::string> >(g_);
	std::vector<std::string> tf = Rcpp::as<std::vector<std::string> >(tf_);
//	gene to tf map
	std::map<std::string,std::vector<std::string> >gtf_map;
	for(int i=0;i<g.size();i++){
		gtf_map[g[i]].push_back(tf[i]);
	}
// now count
	std::map<std::string,std::map<std::string,int> >cmap = get_tf(gtf_map);
	std::map<std::string,double> tfg_idfmap = get_idf(g,tf);
//	return some stuff
	Tfdatamat calc_mat = calc_tfidf(cmap,tfg_idfmap,false,true);
	Rcpp::List out_data = Rcpp::List::create(
			Rcpp::Named("mat")=calc_mat.get_mat(),
			Rcpp::Named("genes")=calc_mat.get_genes(),
			Rcpp::Named("tfs")=calc_mat.get_tfs()
			);
	return Rcpp::wrap(out_data);
}

/*
 * calculate jaccard similarity between the given genes in the gene-tf data.frame
 */
//[[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
SEXP get_ji(SEXP df, int cores){
	Rcpp::DataFrame tf_df(df);
	Rcpp::StringVector g_(tf_df[0]);
	Rcpp::StringVector tf_(tf_df[1]);
	std::vector<std::string> g = Rcpp::as<std::vector<std::string> >(g_);
	std::vector<std::string> tf = Rcpp::as<std::vector<std::string> >(tf_);
//	gene to tf map
	std::map<std::string,std::vector<std::string> >gtf_map;
	for(int i=0;i<g.size();i++){
		gtf_map[g[i]].push_back(tf[i]);
	}
//	get keys (genes) as array
	std::vector<std::string> genes = getkeys(gtf_map);
	int mSize = gtf_map.size();
	Rcpp::NumericMatrix jaccard(mSize,mSize);
//	iterate over map remove duplicates and count overlap
	omp_set_num_threads(cores);
	#pragma omp parallel
	{
		int i,j,c_count,t_count;
		double jaccard_distance;
		std::map<std::string,std::vector<std::string> >::iterator it1,it2;
		std::vector<std::string> tf1,tf2;
		#pragma omp for
		for(i=0;i<mSize;i++){
			tf1 = dedup(gtf_map.find(genes[i])->second);
			jaccard(i,i) = 1.00;
			j=i+1;
			while(j<mSize){
				tf2 = dedup(gtf_map.find(genes[j])->second);
				c_count = getcount(tf1,tf2,true);
				t_count = getcount(tf1,tf2,false);
				jaccard_distance = double(c_count)/double(t_count);
				jaccard(i,j) = jaccard(j,i) = jaccard_distance;
				j++;
			}
		}
	}
	Rcpp::List out_data = Rcpp::List::create(
			Rcpp::Named("mat")=jaccard,
			Rcpp::Named("genes")=genes
			);
	return Rcpp::wrap(out_data);
}
