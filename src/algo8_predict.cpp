// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <math.h>
#include <iostream>
#include <utility>

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>

static double const log2pi = log(2.0 * arma::datum::pi);

using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::export]]
arma::vec replace_cpp(arma::vec x) {

  arma::uword i,nb_clusters;
  arma::vec indices=arma::sort(arma::unique(x));
  nb_clusters=indices.n_elem;

  for(i = 0; i < nb_clusters; i++){
    //arma::vec elements=x.elem(find(x==indices(i)));
    arma::uvec elements=find(x==indices(i));

    arma::vec test=x.elem(find(x==indices(i)));
    arma::uword j,nb_elements_temporaire;
    //nb_elements_temporaire=elements.n_elem;
    nb_elements_temporaire=test.n_elem;
    for(j = 0; j < nb_elements_temporaire; j++){
      //elements(j)=i;
      x(elements(j))=i;
    }
  }
  return(x);
}

int indic(double a, double b){
  int x;
  if(a<b){
    x=1;
  }
  else{
    x=0;
  }
  return(x);
}

arma::vec vector_sub(arma::vec c,arma::uword i) {
  int Ac = c.n_elem;
  arma::vec aRange = arma::regspace<arma::vec>(0, Ac-1);
  arma::vec v = c.elem( find(aRange != i ));
  return v;
}

arma::mat matrix_sub(arma::mat M,arma::uword i,bool row) {
  if(row==TRUE){
    M.shed_row(i);
  }
  else{
    M.shed_col(i);
  }
  return M;
}

// [[Rcpp::export]]
arma::mat dmvt_cpp(arma::mat x,double v,arma::mat mu,arma::mat Sig) {

  arma::uword dim=Sig.n_cols;
  double p=dim;
  //arma::mat densité=lgamma((v+p)/2)-((v+p)/2)*log(1+((x-mu).t()*arma::inv(Sig)*(x-mu))/v)-(lgamma(v/2)+(p/2)*log(v*arma::datum::pi)+0.5*arma::log_det(Sig));
  //arma::mat logdet(1,1);
  double const logdet = sum(log(arma::eig_sym(Sig)));
  arma::mat densité=lgamma((v+p)/2)-((v+p)/2)*log(1+((x-mu).t()*arma::inv(Sig)*(x-mu))/v)-(lgamma(v/2)+(p/2)*log(v*arma::datum::pi)+0.5*logdet);

  return densité;
}

double CRPS_sample(arma::vec sample,double y) {
  arma::uword m=sample.n_elem;

  arma::vec y_ordre=arma::sort(sample);

  double somme=0;
  double m_double=m;

  int i;

  for(i=0;i<m;i++){
    double i_double=i+1;
    double c=0.5;
    somme=somme+((y_ordre(i)-y))*(m_double*(indic(y,y_ordre(i)))-i_double+c);
  }

  double valeur=2*somme/(m*m);


  return valeur;

}

// [[Rcpp::export]]
Rcpp::List partitions_draw_cpp(arma::vec y,arma::mat X_sim,arma::mat X_reg,int M,double alpha,arma::vec c,double v,double s,arma::vec sigma,arma::mat beta,arma::mat m_0,arma::mat Phi,double w,arma::vec time_effects,arma::vec space_effects,arma::vec mixed_effects) {

  //arma::wall_clock timer;

  //timer.tic();

  int t,i;
  arma::uword p_sim = X_sim.n_cols;
  arma::uword p_reg = X_reg.n_cols;
  arma::uword n = X_reg.n_rows;
  arma::vec sigma_nouveau(M,arma::fill::randu);
  arma::mat beta_nouveau(p_reg+1,M);
  arma::mat Id(p_reg+1,p_reg+1,arma::fill::eye);
  arma::vec beta_0(p_reg+1,arma::fill::zeros);

  for (i=0; i<n; i++) {
    arma::vec c_moins_i=vector_sub(c,i);
    arma::vec clusters_uniques=arma::unique(c_moins_i);
    arma::uword k_moins=clusters_uniques.n_elem;

    if (any(c_moins_i==c(i))) {
      sigma_nouveau = 1/(arma::randg(M,arma::distr_param(v/2,2/s)));
      for ( t = 0; t < M; t++) {
        beta_nouveau.col(t)=mvnrnd(beta_0, sigma_nouveau(t)*Id, 1);
      }
    }
    else {
      sigma_nouveau(0) = sigma(c(i));
      beta_nouveau.col(0)=beta.col(c(i));
      sigma=vector_sub(sigma,c(i));
      beta.shed_col(c(i));
      c(i)=k_moins+1;
      c=replace_cpp(c);
      arma::uvec idx=arma::regspace<arma::uvec>(1,M-1);
      sigma_nouveau(idx) = 1/(arma::randg(M-1,arma::distr_param(v/2,2/s)));
      for ( t = 1; t < M; t++) {
        beta_nouveau.col(t)=mvnrnd(beta_0, sigma_nouveau(t)*Id, 1);
      }
    }

    //Calcul des probabilités pour assigner le cluster
    arma::vec proba_log(k_moins+M);
    arma::vec y_i=y.row(i);
    arma::mat X_i(1,p_reg+1,arma::fill::ones);
    arma::uvec indices_X_i=arma::regspace<arma::uvec>(1, p_reg);
    X_i(indices_X_i)=X_reg.row(i);
    c_moins_i=vector_sub(c,i);
    arma::mat X_moins_i(n-1,p_sim);
    X_moins_i=matrix_sub(X_sim,i,TRUE);

    int j;
    for( j = 0; j < k_moins; j++) {
      arma::vec c_cluster=c_moins_i.elem(find(c_moins_i==j));
      double n_cluster_moins_i=(c_cluster).n_elem;
      arma::mat X_cluster_moins_i=(X_moins_i).rows(find(c_moins_i==j));
      arma::mat mean_cluster_moins_i=(arma::mean(X_cluster_moins_i,0)).t();
      arma::mat D_cluster_moins_i(p_sim,p_sim,arma::fill::zeros);
      if(n_cluster_moins_i != 1) {
        D_cluster_moins_i=arma::cov(X_cluster_moins_i,1);
      }
      double v_h=w+n_cluster_moins_i;
      arma::mat mu_h=(m_0+n_cluster_moins_i*mean_cluster_moins_i)/(n_cluster_moins_i+1);
      arma::mat Phi_h=Phi+n_cluster_moins_i*D_cluster_moins_i+(n_cluster_moins_i/(n_cluster_moins_i+1))*(mean_cluster_moins_i*mean_cluster_moins_i.t());
      arma::mat log_w_h=log(n_cluster_moins_i)+dmvt_cpp((X_sim.row(i)).t(),v_h-p_sim+1,mu_h,Phi_h*(n_cluster_moins_i+2)/((v_h-p_sim+1)*(n_cluster_moins_i+1)));
      arma::uword idx_i=i;
      arma::vec mean_y=X_i*beta.col(j)+time_effects(idx_i)+space_effects(idx_i)+mixed_effects(idx_i);
      arma::vec sd(1);
      sd=sqrt(sigma(j));

      if(std::isfinite(as_scalar(log_w_h))==FALSE){
        cout << "eigen: " << arma::eig_sym(Phi_h*(2/(v_h-p_sim+1))) << endl;
        cout << "log_eigen: " << log(arma::eig_sym(Phi_h*(2/(v_h-p_sim+1)))) << endl;
        cout << "D: " << D_cluster_moins_i.is_symmetric() << endl;
        cout << "X'X: " << (mean_cluster_moins_i*mean_cluster_moins_i.t()).is_symmetric() << endl;
        cout << "determinant: " << arma::det(Phi_h*(2/(v_h-p_sim+1))) << endl;
      }
      arma::vec proba_j=arma::log_normpdf(y_i,mean_y,sd)+log_w_h;
      arma::uword idx_y=j;
      proba_log(idx_y)=as_scalar(proba_j);
    }

    arma::mat log_w_0=log(alpha)+dmvt_cpp((X_sim.row(i)).t(),w-p_sim+1,m_0,Phi*(2/(w-p_sim+1)));
    for( j = 0; j < M; j++){
      arma::vec mean_y=X_i*beta_nouveau.col(j)+time_effects(i)+space_effects(i)+mixed_effects(i);
      arma::vec sd(1);
      sd=sqrt(sigma_nouveau(j));
      arma::uword idx_y_2=k_moins+j;
      arma::vec proba_j=arma::log_normpdf(y_i,mean_y,sd)+log_w_0-log(M);
      proba_log(idx_y_2)=as_scalar(proba_j);
    }

    //On tire l'indicatrice du cluster selon les probas
    arma::vec proba_diff=proba_log-arma::max(proba_log);
    double constante_log=arma::max(proba_log)+log(arma::sum(exp(proba_diff)));
    arma::vec ratio=exp(proba_log-constante_log);
    arma::vec span=arma::regspace(0,k_moins+M-1);
    c(i)=as_scalar(Rcpp::RcppArmadillo::sample(span, 1, FALSE, ratio));

    //On garde les parametres associés à au moins une observation
    arma::vec indices_c=arma::sort(arma::unique(c));
    arma::uvec indices_c_2=arma::conv_to<arma::uvec>::from(indices_c);

    arma::mat beta_combine=arma::join_rows(beta,beta_nouveau);
    arma::mat sigma_combine=arma::join_cols(sigma,sigma_nouveau);
    arma::uword nb_clusters_c=indices_c.n_elem;
    beta=beta_combine.cols(indices_c_2);
    sigma=sigma_combine(indices_c_2);

    //On arrange l'index pour que le code soit cohérent
    if(c(i)>k_moins){
      c=replace_cpp(c);
    }
  }
  //double n_time = timer.toc();
  return Rcpp::List::create(Rcpp::Named("beta") = beta,Rcpp::Named("sigma") = sigma,Rcpp::Named("partition") = c );
}

// [[Rcpp::export]]
Rcpp::List partitions_draw_cpp_fast(arma::vec y,arma::mat X_sim,arma::mat X_reg,int M,double alpha,arma::vec c,double v,double s,arma::vec sigma,arma::mat beta,
                                                arma::mat m_0,arma::mat Phi,double w,arma::vec time_effects,arma::vec space_effects,arma::vec mixed_effects,arma::vec sizes_mem,arma::mat sums_mem,arma::cube products_mem) {

  arma::uword p_sim = X_sim.n_cols;
  arma::uword p_reg = X_reg.n_cols;
  arma::uword n = X_reg.n_rows;
  
  int t,i;
  arma::vec sigma_nouveau(M,arma::fill::randu);
  arma::mat beta_nouveau(p_reg+1,M);
  arma::mat Id(p_reg+1,p_reg+1,arma::fill::eye);
  arma::vec beta_0(p_reg+1,arma::fill::zeros);
  bool singleton;

  arma::vec one(1,arma::fill::ones);
  double v_h;
  arma::mat mu_h(1,p_sim);
  arma::mat Phi_h(p_sim,p_sim);
  arma::mat log_w_h(1,1);
  arma::vec sd(1);
  arma::vec mean_y(1);

  for (i=0; i<n; i++) {
    double random_effect=time_effects(i)+space_effects(i)+mixed_effects(i);
    arma::vec c_moins_i=vector_sub(c,i);
    arma::vec clusters_uniques=arma::unique(c_moins_i);
    arma::uword k_moins=clusters_uniques.n_elem;

    //Si pas singleton, tirer M new params. Sinon, deplacer cluster du singleton derniere position et tirer M-1 new params.
    if (any(c_moins_i==c(i))) {
      sigma_nouveau = 1/(arma::randg(M,arma::distr_param(v/2,2/s)));
      for ( t = 0; t < M; t++) {
        beta_nouveau.col(t)=mvnrnd(beta_0, sigma_nouveau(t)*Id, 1);
      }
    }
    else {
      singleton=TRUE;
      sigma_nouveau(0) = sigma(c(i));
      beta_nouveau.col(0)=beta.col(c(i));
      sigma.shed_row(c(i));
      beta.shed_col(c(i));
      sizes_mem.shed_row(c(i));
      sums_mem.shed_col(c(i));
      products_mem.shed_slice(c(i));
      c(i)=k_moins+1;
      c=replace_cpp(c);
      arma::uvec idx=arma::regspace<arma::uvec>(1,M-1);
      sigma_nouveau(idx) = 1/(arma::randg(M-1,arma::distr_param(v/2,2/s)));
      for ( t = 1; t < M; t++) {
        beta_nouveau.col(t)=mvnrnd(beta_0, sigma_nouveau(t)*Id, 1);
      }
    }

    //Calcul des probabilités pour assigner le cluster
    arma::vec proba_log(k_moins+M);
    arma::vec y_i=y.row(i);
    arma::mat X_i(1,p_reg+1,arma::fill::ones);
    arma::uvec indices_X_i=arma::regspace<arma::uvec>(1, p_reg);
    X_i(indices_X_i)=X_reg.row(i);

    int j;
    for( j = 0; j < k_moins; j++) {

      if(j==c(i)){
        int n_modif=sizes_mem(j)-1;
        arma::mat sum_modif=sums_mem.col(j)-(X_sim.row(i)).t();
        arma::mat sum_carre_modif=products_mem.slice(j)-(X_sim.row(i).t())*X_sim.row(i);
        arma::mat mean_modif=(sum_modif/n_modif);
        arma::mat D_modif=sum_carre_modif-(sum_modif*(sum_modif).t())/(n_modif);
        v_h=w+n_modif;
        mu_h=(m_0+n_modif*mean_modif)/(n_modif+1);
        Phi_h=Phi+D_modif+(n_modif/(n_modif+1))*((mean_modif-m_0)*(mean_modif-m_0).t());
        //log_w_h(j)=as_scalar(log(n_modif)+dmvt_cpp((X_sim.row(i)).t(),v_h-p_sim+1,mu_h,Phi_h*(n_modif+2)/((v_h-p_sim+1)*(n_modif+1))));
        sd=sqrt(sigma(j));
        mean_y=X_i*beta.col(j)+random_effect;
        proba_log(j)=as_scalar(log(n_modif)+dmvt_cpp((X_sim.row(i)).t(),v_h-p_sim+1,mu_h,Phi_h*(n_modif+2)/((v_h-p_sim+1)*(n_modif+1)))+arma::log_normpdf(y_i,mean_y,sd));

      } else{
        int n_h=sizes_mem(j);
        arma::mat mean_h=sums_mem.col(j)/n_h;
        arma::mat D_h=products_mem.slice(j)-(sums_mem.col(j)*(sums_mem.col(j)).t())/(n_h);
        v_h=w+n_h;
        mu_h=(m_0+n_h*mean_h)/(n_h+1);
        Phi_h=Phi+D_h+(n_h/(n_h+1))*((mean_h-m_0)*(mean_h-m_0).t());
        sd=sqrt(sigma(j));
        mean_y=X_i*beta.col(j)+random_effect;
        proba_log(j)=as_scalar(log(n_h)+dmvt_cpp((X_sim.row(i)).t(),v_h-p_sim+1,mu_h,Phi_h*(n_h+2)/((v_h-p_sim+1)*(n_h+1)))+arma::log_normpdf(y_i,mean_y,sd));;
      }
    }

    arma::mat log_w_0=log(alpha)+dmvt_cpp((X_sim.row(i)).t(),w-p_sim+1,m_0,Phi*(2/(w-p_sim+1)));
    for( j = 0; j < M; j++){
      arma::vec mean_y=X_i*beta_nouveau.col(j)+random_effect;
      arma::vec sd(1);
      sd=sqrt(sigma_nouveau(j));
      arma::uword idx_y_2=k_moins+j;
      arma::vec proba_j=arma::log_normpdf(y_i,mean_y,sd)+log_w_0-log(M);
      proba_log(idx_y_2)=as_scalar(proba_j);
    }

    //On tire l'indicatrice du cluster selon les probas
    arma::vec proba_diff=proba_log-arma::max(proba_log);
    double constante_log=arma::max(proba_log)+log(arma::sum(exp(proba_diff)));
    arma::vec ratio=exp(proba_log-constante_log);
    arma::vec span=arma::regspace(0,k_moins+M-1);
    int c_old=c(i);
    c(i)=as_scalar(Rcpp::RcppArmadillo::sample(span, 1, FALSE, ratio));

    //On garde les parametres associés à au moins une observation
    arma::vec indices_c=arma::sort(arma::unique(c));
    arma::uvec indices_c_2=arma::conv_to<arma::uvec>::from(indices_c);
    arma::mat beta_combine=arma::join_rows(beta,beta_nouveau);
    arma::mat sigma_combine=arma::join_cols(sigma,sigma_nouveau);
    beta=beta_combine.cols(indices_c_2);
    sigma=sigma_combine(indices_c_2);

    //On arrange l'index pour que le code soit cohérent
    if(c(i)>k_moins){
      c=replace_cpp(c);
    }

    //Si s_i changed, mettre a jour les quantite dans cluster de depart et d'arrive en ignorant mise a jour si cluster de départ si singleton.
    //On enleve la contribution dans l'ancien cluster. On ajoute la contribution dans le nouveau cluster. SI c'est un brand new cluster, on ajoute quantité simple , sinon on modifie anciennes quantités.
    if(singleton!=TRUE){
      sizes_mem.row(c_old)=sizes_mem(c_old)-1;
      sums_mem.col(c_old)=sums_mem.col(c_old)-(X_sim.row(i)).t();
      products_mem.slice(c_old)=products_mem.slice(c_old)-(X_sim.row(i).t())*X_sim.row(i);
    }
    if(c(i)>=k_moins){
      sizes_mem.insert_rows(sizes_mem.n_rows,one);
      sums_mem.insert_cols(sums_mem.n_cols,(X_sim.row(i)).t());
      products_mem.insert_slices(products_mem.n_slices,(X_sim.row(i).t())*X_sim.row(i));
    } else{
      sizes_mem.row(c(i))=sizes_mem.row(c(i))+1;
      sums_mem.col(c(i))=sums_mem.col(c(i))+(X_sim.row(i)).t();
      products_mem.slice(c(i))=products_mem.slice(c(i))+(X_sim.row(i).t())*X_sim.row(i);
    }
  }
  return Rcpp::List::create(Rcpp::Named("beta") = beta,Rcpp::Named("sigma") = sigma,Rcpp::Named("partition") = c,Rcpp::Named("sizes_mem") = sizes_mem,Rcpp::Named("sums_mem") = sums_mem,Rcpp::Named("products_mem") = products_mem );
}

// [[Rcpp::export]]
Rcpp::List predict_cpp_new(int obs,int station,int time,arma::mat X_train_sim,arma::mat x_futur_sim,arma::mat x_futur_reg,double y_test,Rcpp::List beta_list ,Rcpp::List sigma_list,Rcpp::List space_effects_list,Rcpp::List time_effects_list,Rcpp::List mixed_effects_list,
                                        Rcpp::List clusters_sizes_list,Rcpp::List partition_list,arma::vec alpha_vec,arma::uword nb_density,int nb_draw,double w,arma::mat Phi,double v,double s, arma::mat m_0) {


  arma::uword p_sim = X_train_sim.n_cols;
  arma::uword p_reg = x_futur_reg.n_elem;
  arma::vec beta_0(p_reg+1,arma::fill::zeros);

  //On ajoute constante 1 au début du vecteur x_futur
  arma::mat x_futur_1(p_reg + 1, 1, arma::fill::ones);
  arma::uvec indices=arma::regspace<arma::uvec>(1, p_reg);
  x_futur_1(indices)=x_futur_reg;
  arma::mat beta_densite(p_reg + 1,nb_density);
  arma::vec sigma_densite(nb_density);
  arma::vec mean_futur(nb_density);
  arma::vec mu_h(p_sim);
  arma::mat Phi_h(p_sim,p_sim);

  arma::uword b;
  for(b=0; b < nb_density; b++){

    arma::vec taille_clusters_b=clusters_sizes_list[b];
    arma::vec partition_b=partition_list[b];
    arma::vec effet_spatial_b=space_effects_list[b];
    arma::vec effet_temps_b=time_effects_list[b];
    arma::vec effet_mixed_b=mixed_effects_list[b];
    arma::mat beta_b=beta_list[b];
    arma::vec sigma_b=sigma_list[b];
    double alpha_b=alpha_vec(b);
    arma::uword nb_clusters=taille_clusters_b.n_elem;
    arma::vec w_h_futur(nb_clusters+1);

    //Calcul du poids pour que observation soit dans nouveau cluster
    w_h_futur(nb_clusters)=arma::as_scalar(log(alpha_b)+dmvt_cpp(x_futur_sim,w-p_sim+1,m_0,2*Phi/(w-p_sim+1)));
    arma::uword h;
    arma::vec mu_h(p_sim);
    arma::mat Phi_h(p_sim,p_sim);

    //Calcul des poids pour que observations soit dans cluster existant
    for(h=0; h<nb_clusters; h++){
      int n_h=taille_clusters_b(h);
      arma::mat X_h(n_h,p_sim);
      arma::uvec cluster_h=arma::find(partition_b==h);
      X_h=X_train_sim.rows(cluster_h);
      arma::vec mean_h(p_sim);
      mean_h=(arma::mean(X_h)).t();

      int q;
      arma::mat D_h(p_sim,p_sim);
      for (q=0; q<n_h; q++) {
        D_h=D_h+(((X_h.row(q)).t()-mean_h)*(((X_h.row(q)).t()-mean_h)).t())/n_h;
      }
      double v_h=w+n_h;
      mu_h=(m_0+n_h*mean_h)/(n_h+1);
      Phi_h=Phi+n_h*D_h+(n_h/(n_h+1))*((mean_h-m_0)*(mean_h-m_0).t());
      w_h_futur(h)=log(n_h)+arma::as_scalar(dmvt_cpp(x_futur_sim,v_h-p_sim+1,mu_h,2*Phi_h/(v_h-p_sim+1)));

    }
    double w_max=w_h_futur.max();
    double constante_w=w_max+log(arma::sum(exp(w_h_futur-w_max)));
    arma::vec ratio(nb_clusters+1);
    ratio=exp(w_h_futur-constante_w);
    arma::vec span=arma::regspace(0,nb_clusters);
    int s_futur=arma::as_scalar(Rcpp::RcppArmadillo::sample(span, 1, FALSE, ratio));
    arma::vec beta_futur(p_reg+1);
    arma::mat Id(p_reg+1,p_reg+1,arma::fill::eye);

    if(s_futur==(nb_clusters)){
      double sigma_futur=1/(arma::as_scalar(arma::randg(1,arma::distr_param(v/2,2/s))));
      beta_futur=mvnrnd(beta_0,sigma_futur*Id,1);
      beta_densite.col(b)=beta_futur;
      sigma_densite(b)=sigma_futur;
      std::cout << "nouveau cluster pour observation " << obs+1 << endl;
    } else{
      double sigma_futur=sigma_b(s_futur);
      beta_futur=(beta_b).col(s_futur);
      beta_densite.col(b)=beta_futur;
      sigma_densite(b)=sigma_futur;
    }
    mean_futur(b)=arma::as_scalar(((x_futur_1).t())*(beta_futur)+effet_spatial_b(station)+effet_temps_b(time)+effet_mixed_b(obs));
  }

  int j;
  arma::vec sample_predi(nb_draw);
  for(j=0;j<nb_draw;j++){
    arma::vec span_tir=arma::regspace(0,nb_density-1);
    int tir_b=arma::as_scalar(Rcpp::RcppArmadillo::sample(span_tir, 1, FALSE));
    arma::mat beta_tir=beta_densite.col(tir_b);
    double sigma_tir=sigma_densite(tir_b);
    arma::vec effet_spatial_tir=space_effects_list[tir_b];
    arma::vec effet_temps_tir=time_effects_list[tir_b];
    arma::vec effet_mixed_tir=mixed_effects_list[tir_b];
    double mean_tir=arma::as_scalar(((x_futur_1).t())*(beta_tir)+effet_spatial_tir(station)+effet_temps_tir(time)+effet_mixed_tir(obs) );
    sample_predi(j)=arma::as_scalar(arma::randn(arma::distr_param( mean_tir, sqrt(sigma_tir) ) ) );
  }
  if((obs+1)%100==0){
    cout << "predicting observation: " << obs+1 << endl;
  }
  return Rcpp::List::create(Rcpp::Named("predictions") = arma::mean(mean_futur),Rcpp::Named("CRPS") = CRPS_sample(sample_predi,y_test),Rcpp::Named("sample") =sample_predi);
}
