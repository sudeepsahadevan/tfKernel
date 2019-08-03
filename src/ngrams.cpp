#include <Rcpp.h>
#include <omp.h>

// get n grams for real and remove duplicates
std::vector<std::string> get_ngrams(std::string instr,int ng,bool dedup){
	std::vector<std::string> ngram_vec;
	int len = instr.length();
	int i = 0;
	while(i+ng<=len){
		std::string s = instr.substr(i,ng);
		ngram_vec.push_back(s);
		i++;
	}
// remove duplicate grams
// found here: http://stackoverflow.com/questions/1041620/whats-the-most-efficient-way-to-erase-duplicates-and-sort-a-vector
	if(dedup){
		std::set<std::string> set_n;
		unsigned size = ngram_vec.size();
		for( unsigned i = 0; i < size; ++i ) set_n.insert( ngram_vec[i] );
		ngram_vec.assign( set_n.begin(), set_n.end() );
	}
	return ngram_vec;
}
// find common stuff
std::vector<std::string> get_common(std::vector<std::string> n_g1,std::vector<std::string>n_g2){
	std::sort(n_g1.begin(),n_g1.end());
	std::sort(n_g2.begin(),n_g2.end());
	std::vector<std::string> common;
	std::set_intersection(n_g1.begin(),n_g1.end(),n_g2.begin(),n_g2.end(),std::back_inserter(common));
	return common;
}

class NgramData{
private:
	std::string id1;
	std::string id2;
	std::vector<std::string> cgrams;
	int clen;
public:
	NgramData();
	NgramData(std::string id1,std::string id2,std::vector<std::string> cgrams,int clen);
	NgramData(std::string id1,std::string id2,int clen);
//	 setters
	/**
	void set_id1();
	void set_id2();
	void set_cgrams();
	void set_clen();
	**/
//	getters
	std::string get_id1();
	std::string get_id2();
	std::vector<std::string> get_cgrams();
	int get_clen();
};

NgramData::NgramData(std::string id1,std::string id2,std::vector<std::string> cgrams,int clen){
	this->id1=id1;
	this->id2=id2;
	this->cgrams=cgrams;
	this->clen=clen;
}

NgramData::NgramData(std::string id1,std::string id2,int clen){
	this->id1=id1;
	this->id2=id2;
	this->clen=clen;
}

std::string NgramData::get_id1(){
	return id1;
}

std::string NgramData::get_id2(){
	return id2;
}

std::vector<std::string> NgramData::get_cgrams(){
	return cgrams;
}

int NgramData::get_clen(){
	return clen;
}

// Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
// Sys.setenv("PKG_LIBS"="-fopenmp")

// [[Rcpp::export]]
SEXP ngrams(SEXP strs,SEXP n,SEXP rmdup){
BEGIN_RCPP
	std::string instr = Rcpp::as<std::string>(strs);
	int ng = Rcpp::as<int>(n);
	bool dup_rm = Rcpp::as<bool>(rmdup);
	return Rcpp::wrap(get_ngrams(instr,ng,dup_rm));
END_RCPP
}

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
SEXP ngrams_df(SEXP df,SEXP n,SEXP rmdup,SEXP cores){
BEGIN_RCPP
	Rcpp::DataFrame seq_df = Rcpp::DataFrame(df);
	int ng = Rcpp::as<int>(n);
	bool dup_rm = Rcpp::as<bool>(rmdup);
	int threads = Rcpp::as<int>(cores);
	// do something to the data.frame
	Rcpp::StringVector ath_id = seq_df[0];
	Rcpp::StringVector seq = seq_df[1];
	std::vector<std::string> id_vec = Rcpp::as< std::vector< std::string> >(ath_id);
	std::vector<std::string> seq_vec = Rcpp::as< std::vector<std::string> >(seq);
	// put them into an array
	std::vector<std::vector<std::string> >  ngram_comm;
	for(int i=0;i<seq_vec.size();i++){
		std::vector<std::string> ngram_tokens = get_ngrams(seq_vec[i],ng,dup_rm);
		ngram_comm.push_back(ngram_tokens);
	}
	std::vector<NgramData> ng_vector;
	int k,m;
	int in_a=0;
	std::vector<std::string> k_gram,m_gram,common_gram;
	std::string id1,id2;
	std::vector<NgramData> innerVec;
	omp_set_num_threads(threads);
	#pragma omp parallel shared(ngram_comm,ng_vector)
	{
		#pragma omp for private(k,k_gram,id1,m,m_gram,id2,common_gram,innerVec) schedule(dynamic)
		for(k=0;k<ngram_comm.size();k++){
			k_gram = ngram_comm[k];
			id1 = id_vec[k];
			for(m=0;m<ngram_comm.size();m++){
				if(m>=k){
					m_gram = ngram_comm[m];
					id2 = id_vec[m];
					common_gram = get_common(k_gram,m_gram);
					if(common_gram.size()>0){
						NgramData ng_d(id1,id2,common_gram.size());
						innerVec.push_back(ng_d);
					}
				}
			}
			#pragma omp critical
			{
				ng_vector.insert(ng_vector.end(),innerVec.begin(),innerVec.end());
				innerVec.clear();
			}
		}
	}
	// vectors to generate data.frame()
	std::vector<std::string> id1_vec;
	std::vector<std::string> id2_vec;
	std::vector<int> len_vec;
	for(int n=0;n<ng_vector.size();n++){
		NgramData ng = ng_vector[n];
		id1_vec.push_back(ng.get_id1());
		id2_vec.push_back(ng.get_id2());
		len_vec.push_back(ng.get_clen());
	}
	Rcpp::DataFrame retFrame = Rcpp::DataFrame::create(
			Rcpp::Named("id1")=id1_vec,
			Rcpp::Named("id2")=id2_vec,
			Rcpp::Named("ngrams")=len_vec);
	return Rcpp::wrap(retFrame);
END_RCPP
}


/*

SEXP ngrams_df(SEXP df,SEXP n,SEXP rmdup){
BEGIN_RCPP
Rcpp::DataFrame seq_df = Rcpp::DataFrame(df);
	int ng = Rcpp::as<int>(n);
	bool dup_rm = Rcpp::as<bool>(rmdup);
	int threads = Rcpp::as<int>(cores);
//	 do something to the data.frame
	std::vector<std::string> id_vec = Rcpp::as< std::vector< std::string> >(seq_df[0]);
	std::vector<std::string> seq_vec = Rcpp::as< std::vector<std::string> >(seq_df[1]);
//	 put them into an array
	std::vector<std::vector<std::string> >  ngram_comm;
	for(int i=0;i<seq_vec.size();i++){
		std::vector<std::string> ngram_tokens = get_ngrams(seq_vec[i],ng,dup_rm);
		ngram_comm.push_back(ngram_tokens);
	}
//	vector for storing return data
	std::vector<NgramData> ng_vector;
	for(int k =0;k<ng_vector.size();k++){
		k_gram = ngram_comm[k];
		id1 = id_vec[k];
		int m=k;
		while(m<ng_vector.size()){
			m_gram = ngram_comm[m];
			id2 = id_vec[m];
			common_gram = get_common(k_gram,m_gram);
			if(common_gram.size()>0){
				NgramData ng_d(id1,id2,common_gram.size());
				ng_vector.push_back(ng_d);
			}
			m++;
		}
	}
//	 vectors to generate data.frame()
	std::vector<std::string> id1_vec;
	std::vector<std::string> id2_vec;
	std::vector<int> len_vec;
	for(int n=0;n<ng_vector.size();n++){
		NgramData ng = ng_vector[n];
		id1_vec.push_back(ng.get_id1());
		id2_vec.push_back(ng.get_id2());
		len_vec.push_back(ng.get_clen());
		}
	Rcpp::DataFrame dataframe = Rcpp::DataFrame::create(
			Rcpp::Named("id1")=id1_vec,
			Rcpp::Named("id2")=id2_vec,
			Rcpp::Named("ngrams")=len_vec);
	return Rcpp::wrap(dataframe);
END_RCPP
}
*/
