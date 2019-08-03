#include <Rcpp.h>
#include <omp.h>
#include <typeinfo>
/*
 * Use openmpi for fast matrix multiplication
 */

// Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
// Sys.setenv("PKG_LIBS"="-fopenmp")

/*
 * do matrix multiplication using n cores
 */

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
SEXP mat_mult(SEXP m1_,SEXP m2_,int cores){
	Rcpp::NumericMatrix m1(m1_);
	Rcpp::NumericMatrix m2(m2_);
//	a few ints
	int mrow = m1.nrow();
	int mcol = m1.ncol();
	int ncol = m2.ncol();
//	new matrix
	Rcpp::NumericMatrix mprod (mrow, ncol);
	omp_set_num_threads(cores);
	#pragma omp parallel
	{
		int i,j,k;
		#pragma omp for
		for(i=0;i<mrow;i++){
			for(j=0;j<ncol;j++){
				for(k=0;k<mcol;k++){
					mprod(i,j)+=m1(i,k)*m2(k,j);
				}
			}
		}
	}
	return(mprod);
}


/*
 * estimate sigma/gamma for rbf/laplacian kernels using n cores
 */

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
std::vector<double> est_sig(SEXP mat_,int cores){
	Rcpp::NumericMatrix mat(mat_);
	int mrow = mat.nrow();
	int mcol = mat.ncol();
//	define a vector for distances
	std::vector<double> dist;
	omp_set_num_threads(cores);
	#pragma omp parallel
	{
		int i,j,k;
		double sum;
		std::vector<double> dist_private;
		#pragma omp for
		for(i = 0;i<mrow;i++){
			j = i+1;
			while(j<mrow){
				sum = 0;
				for(k=0;k<mcol;k++){
					sum+= pow(double(mat(i,k)-mat(j,k)),2);
				}
				dist_private.push_back(sum);
				j++;
			}
		}
		#pragma omp critical
		dist.insert(dist.end(),dist_private.begin(),dist_private.end());
	}
	return(dist);
}


/*
 * calculate laplacian/rbf kernel using n cores
 */

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
SEXP rbf_kern(SEXP mat_,double sigma,int cores,bool lap=false){
	Rcpp::NumericMatrix mat(mat_);
	int mrow = mat.nrow();
	int mcol = mat.ncol();
//	define a new matrix to return
	Rcpp::NumericMatrix lk(mrow,mrow);
	omp_set_num_threads(cores);
	#pragma omp parallel
	{
		int i,j,k;
		double sum,lk_e;
		#pragma omp for
		for(i=0;i<mrow;i++){
			j=i+1;
			lk(i,i)=1.00;
			while(j<mrow){
				sum = 0;
				for(k=0;k<mcol;k++){
					sum+=double(pow((mat(i,k)-mat(j,k)),2));
				}
				if(lap){
					sum = sqrt(sum);
				}
				lk_e = double(exp(-sigma*sum));
				lk(i,j)=lk(j,i)=lk_e;
				j++;
			}
		}
	}
	return(lk);
}

/*
 * for a given matrix and vector of distance values, calculate rbf/laplacian kernel
 * @param mat_ : numeric matrix
 * @param vec_ : distance vector for each data point
 * @param cores : number of cores to use
 * @param lap : boolean, whether to calculate laplacian or rbf kernel
 */

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
SEXP rbf_kern_nn(SEXP mat_,SEXP vec_,int cores = 3, bool lap = false){
	Rcpp::NumericMatrix mat(mat_);
	int mrow = mat.nrow();
	int mcol = mat.ncol();
	Rcpp::NumericVector vec(vec_);
	std::vector<double> dv = Rcpp::as<std::vector<double> >(vec);
//	matrix to return
	Rcpp::NumericMatrix nmat(mrow,mrow);
	omp_set_num_threads(cores);
	#pragma omp parallel default(none) shared(dv,mat,mrow,mcol,lap,nmat)
	{
		int i,j,k;
		double sum,dv_i,dv_j,dv_sum;
		#pragma omp for
		for(i=0;i<mrow;i++){
			nmat(i,i)=1.00;
			j = i;
			dv_i = dv[i];
			while(j<mrow){
				sum = 0.0;
				dv_j = dv[j];
				dv_sum = double((dv_i+dv_j)/2);
				nmat(i,i) = 1.00;
				for(k=0;k<mcol;k++){
					sum += pow(double(mat(i,k)-mat(j,k)),2);
				}
				if(lap){
					sum = sqrt(sum)/sqrt(dv_sum);
				}else{
					sum = sum/(dv_sum);
				}
				sum = exp(-sum);
				nmat(i,j) = nmat(j,i) = sum;
				j++;
			}
		}
	}
	return(nmat);
}

/*
 * calculate euclidian distance
 * given an input matrix calculate the euclidean distance
 * @param mat_ : numeric matrix
 * @param cores: number of cores to use
 * @param squared: boolean, to return the squared distance or not
 */

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
SEXP dist_mat(SEXP mat_,int cores, bool squared = true){
	Rcpp::NumericMatrix mat(mat_);
	int mrow = mat.nrow();
	int mcol = mat.ncol();
	Rcpp::NumericMatrix dimat(mrow,mrow);
	omp_set_num_threads(cores);
	#pragma omp parallel
	{
		int i,j,k;
		double sum;
		for(i=0;i<mrow;i++){
			dimat(i,i)=0.0;
			j = i;
			while(j<mrow){
				sum = 0.0;
				for(k=0;k<mcol;k++){
					sum += pow(double(mat(i,k)-mat(j,k)),2);
				}
				if(!squared){
					sum = sqrt(sum);
				}
				dimat(i,j)=dimat(j,i)=sum;
				j++;
			}
		}
	}
//	std::cout<<typeid(dimat).name()<<"\n";
	return(dimat);
}

/*
 * given a distance matrix, sort the matrix and get distances for n nearest neighbors
 * @param dist: euclidean distance matrix
 * @param nn: N nearest neighbors
 * // [[Rcpp::export]]
 */

SEXP get_nn_dist(SEXP dist_,int nn = 10){
	Rcpp::NumericMatrix dist(dist_);
	int mrow = dist.nrow();
	int mcol = dist.ncol();
	Rcpp::NumericMatrix distnn(mrow,nn);
	std::vector<double>::iterator it;
	for(int i=0;i<mrow;i++){
		Rcpp::NumericVector v(dist.row(i));
		std::vector<double> dv = Rcpp::as<std::vector<double> >(v);
		std::sort(dv.begin(),dv.end());
		std::vector<double> vnn (dv.begin()+1,dv.begin()+(nn+1));
		int j = 0;
		for(it = vnn.begin();it!=vnn.end();it++){
			distnn(i,j) = *it;
			j++;
		}
	}
	return(distnn);
}

/*
 * given a matrix, calculate the  euclidean|squared distance and
 * return the top n nearest neighbor values
 * @param mat_ : numeric matrix
 * @param cores: number of cores to use
 * @param squared: boolean, to return the squared distance or not
 * @param dist: euclidean distance matrix
 * @param nn: N nearest neighbors
 */
// [[Rcpp::export]]
SEXP calc_nn(SEXP mat_,int cores = 3, bool squared = false, int nn = 10){
	Rcpp::NumericMatrix dm(dist_mat(mat_,cores,squared));
	return(get_nn_dist(dm,nn));
}

/*
 * calculate linear kernel
 */

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
SEXP linear_kern(SEXP mat_,int cores){
	Rcpp::NumericMatrix mat(mat_);
	int mrow = mat.nrow();
	int mcol = mat.ncol();
//	define a new matrix to return
	Rcpp::NumericMatrix link(mrow,mrow);
	omp_set_num_threads(cores);
	#pragma omp parallel
	{
		int i,j,k;
		double sumab;
		for(i=0;i<mrow;i++){
			j = i;
			while(j<mrow){
				sumab = 0.0;
				for(k=0;k<mcol;k++){
					sumab+=mat(i,k)*mat(j,k);
				}
				link(i,j)=link(j,i)=sumab;
				j++;
			}
		}
	}
	return(link);
}

/*
 * calculate cosine similarity for a given matrix using n cores
 */

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
SEXP cosine_kern(SEXP mat_,int cores){
	Rcpp::NumericMatrix mat(mat_);
	int mrow = mat.nrow();
	int mcol = mat.ncol();
//	define a new matrix to return
	Rcpp::NumericMatrix cosk(mrow,mrow);
	omp_set_num_threads(cores);
	#pragma omp parallel
	{
		int i,j,k;
		double sumab,sumaa,sumbb;
		#pragma omp for
		for(i=0;i<mrow;i++){
			j = i+1;
			cosk(i,i)=1.00;
			while(j<mrow){
				sumab = sumaa = sumbb = 0;
				for(k=0;k<mcol;k++){
					sumab += mat(i,k)*mat(j,k);
					sumaa += pow(mat(i,k),2);
					sumbb += pow(mat(j,k),2);
				}
				cosk(i,j)=cosk(j,i)=sumab/(sqrt(sumaa)*sqrt(sumbb));
				j++;
			}
		}
	}
	return(cosk);
}

