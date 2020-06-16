#include <cmath>
#include <cstdio>
#include <vector>
#include <time.h>
#include <sys/stat.h>
using vd = std::vector<double>;
using vvd = std::vector<vd>;
using vvvd = std::vector<vvd>;
using vvvvd = std::vector<vvvd>;


/**********************config*************************/
// terminal に出力して時間確認 
const int OUTPUT_TERMINAL = 1;
const int DEBUG_f = 0;
// 初期条件 -> 0:長方形, 1:sinx 
const int INITIAL = 1;
// mesh -> 0:linear, 1:log
const int LOGMESH = 0;
/*****************************************************/

/******************計算条件********************/
const double Lmax = 10.0;
const double Lmin = 1e-1;
const double Lr = 1.0;
const int NR = 100 + 1;
const double dr = Lmax/(NR-1);
const double kappa = 1.0;
const double OMEGA = 1.8;
const double dt = 0.001;
const double TMAX = 1.0 + 1e-9;
const double EPS = 1e-10;

vd r(NR),rp(NR),rm(NR);
void grid_set();

/*********************************************/
/*
以下出力枚数の調整用
時刻 0 ~ endtime の間のプロファイルを時間 DT ごとに出力させるための変数達
TIME に配列 (DT, 2DT, ..., ENDTIME) をsetする
*/
const double T_EPS = 1.0e-10;
const double DT = 0.01;
const double ENDTIME = 1.0;
vd TIME;
//出力時刻を set する関数
void TIME_set();
/*********************************************/


/*************************** diffusion を解く関数 ***************************/
// void init(vvd &f);
void init(vd &f);
void boundary(vd &nf);
int diffusion(vd &f);
double SOR(vd &f, vd &a);
/*************************** diffusion を解く関数 ***************************/

/*********************************ファイル関連*********************************/
FILE *condition_fp = fopen("../data/condition.csv", "w");
FILE *time_fp = fopen("../data/output_time.csv", "w");
FILE *conservation_fp = fopen("../data/conservation.csv", "w");
FILE *f_fp;
char f_filename[50];
// その時刻における f の値をファイルにアウトプット
void output(double t, vd &f);
void conservation(double t, vd &f);
/*****************************************************************************/

int main(){
  // 時間計測
  clock_t start_t, end_t;
  start_t = time(NULL);

  vd f(NR+1);
  int ti = 0;
  double nf;
  double t = 0.0;
  TIME_set();
  grid_set();
  // for(int i = 0; i < NR; i++){
  //   if(i < NR-1) printf("%f ", r[i]);
  //   else printf("%f\n", r[i]);
  // }
  // return 0;
  init(f);
  fprintf(time_fp, "time\n");
  fprintf(conservation_fp, "time,F\n");
  fprintf(condition_fp, "Lmax,Lr,dt,dr,kappa,NR\n%e,%e,%e,%e,%e,%d\n",Lmax,Lr,dt,dr,kappa,NR);
  output(t, f);
  conservation(t, f);
  printf("****************CALICULATION START****************\n");
  while(t < TMAX) {
    // f から nf を求める
    if(diffusion(f) == -1){
      printf("Diffusion equation does not converge!!!\n");
      return 0;
    }
    boundary(f);
    t += dt;
    // nf が求め終わったのでこれを出力
    if(ti < TIME.size() && t > TIME[ti] - T_EPS){
      output(t, f);
      conservation(t, f);
      ti++;
    }
  }

  fclose(time_fp);
  fclose(condition_fp);
  fclose(conservation_fp);
  printf("******************NORMAL END******************\n");
  end_t = time(NULL);
  printf("This calculatioin took %ld second \n", end_t - start_t);
  return 0;
}

void init(vd &f){
  if(INITIAL == 0){
    for(int i = 0; i < NR; i++) {
      if(r[i] < Lr + EPS){
        f[i] = 1.0;
      }
      else{
        f[i] = 0.0;
      }
    }
  }

  if(INITIAL == 1){
    for(int i = 0; i < NR; i++) {
      f[i] = (Lr/M_PI)*std::sin(M_PI*r[i]/Lr)/r[i];
    }
  }
  // boundary(f);
}

void TIME_set(){
  double t = DT;
  while(t < ENDTIME + T_EPS){
    TIME.push_back(t);
    t += DT;
  }
}

void grid_set(){
  if(LOGMESH){
    double logdr = (std::log10(Lmax)-std::log10(Lmin))/(NR-1);
    for(int i = 0; i < NR; i++) {
      r[i] = Lmin*std::pow(10.0,i*logdr);
      
    }
    for(int i = 0; i < NR-1; i++) {
      rp[i] = std::sqrt(r[i]*r[i+1]);
    }
    for(int i = 1; i < NR; i++) {
      rm[i] = rp[i-1];
    }
    rm[0] = 0.0;
  }
  else{
    for(int i = 0; i < NR; i++) {
      r[i] = dr*i;
    }
    for(int i = 0; i < NR-1; i++) {
      rp[i] = 0.5*(r[i]+r[i+1]);
    }
    for(int i = 1; i < NR; i++) {
      rm[i] = rp[i-1];
    }
    rm[0] = 0.0;

    // 初期条件設定時の 0 除算を避けるための処方箋
    r[0] = 1e-20;
  }
}


void output(double t, vd &f){
  sprintf(f_filename, "../data/f/%.3f.csv", t);
  FILE *fp = fopen(f_filename, "w");
  fprintf(time_fp, "%f\n", t);

  fprintf(fp, "r,f\n");
  if(OUTPUT_TERMINAL) printf("time:%f\n", t);
  
  for(int i = 0; i < NR; i++) {
    if(DEBUG_f){
      if(i < NR-1) printf("i:%d r:%e f:%f ", i, r[i], f[i]);
      else printf("i:%d r:%e f:%f\n", i, r[i], f[i]);
    }
    fprintf(fp, "%e,%e\n", r[i], f[i]);
  }
  fclose(fp);
}

void boundary(vd &f){
  // こいつの役割が微妙
  f[NR-1] = 0.0;
}

int diffusion(vd &f){
  // f から nf を求める
  vd a(NR);
  int imax = 99999;
  double df;
  for(int i = 0; i < NR-1; i++) {
    a[i] = f[i];
  }
  // nf を求める
  for(int icnt = 0; icnt < imax; icnt++) {
    // f から 仮のnuを求める
    // |f-nd| < EPS ならこの f をアクセプト
    df = SOR(f, a);
    if(df < EPS) return icnt;
  }

  // imax 回反復しても収束しないなら -1 を返して終了
  return -1;
}

double SOR(vd &f, vd &a){
  double nf, df = 0.0, tmp, l, co_p, co_m;
  int i = 0;
  // 仮の f から nf (仮のnu)を求める
  {// i = 0
    l = 3.0*kappa*dt/(rp[i]*rp[i]*rp[i]);
    co_p = rp[i]*rp[i]/(r[i+1]-r[i]);
    tmp = 1.0+l*co_p;
    nf = OMEGA*( l*(co_p*f[i+1]) + a[i] )/tmp;
    nf += (1.0-OMEGA)*f[i];
    df = std::max(std::abs(nf-f[i]), df);
    f[i] = nf;
  }
  for(i = 1; i < NR-1; i++) {
    l = 3.0*kappa*dt/(rp[i]*rp[i]*rp[i]-rm[i]*rm[i]*rm[i]);
    co_p = rp[i]*rp[i]/(r[i+1]-r[i]);
    co_m = rm[i]*rm[i]/(r[i]-r[i-1]);
    tmp = 1.0+l*(co_p + co_m);
    nf = OMEGA*( l*(co_p*f[i+1]+co_m*f[i-1]) + a[i] )/tmp;
    nf += (1.0-OMEGA)*f[i];
    df = std::max(std::abs(nf-f[i]), df);
    f[i] = nf;
  }
  return df;
}

void conservation(double t, vd &f){
  double sum = 0.0;
  for(int i = 0; i < NR; i++) {
    double vol = 4.0*M_PI/3.0*(rp[i]*rp[i]*rp[i]-rm[i]*rm[i]*rm[i]);
    sum += f[i]*vol;
  }
  if(OUTPUT_TERMINAL) printf("conservation:%e\n", sum);
  fprintf(conservation_fp, "%e,%e\n",t,sum);
}
