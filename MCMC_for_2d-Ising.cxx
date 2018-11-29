# include <iostream>
# include <math.h>
# include <iomanip>
# include <random>
# include <fstream>
# include <time.h>
# include <stdio.h>
# include <limits.h>

//This calulate the specific_heat of the 2-d square Ising model whose lattice-size is 50*50.

using namespace std;

#define N 50 //system size
#define TIME 100000000
#define step 6
#define seed 3
#define bin 10
#define T_i 2.25
#define T_f 2.32
#define burn_in 2000

#define J 1.0//coupling constant

int a[N][N];

void initial();
void monte_iti_even(double &,double &,double &);
void monte_iti_even_non_energy(double &,double &);
void energy_cal(double &);
void work_even(double &,double &,double &);
void rand_set();
void rand_reset();


mt19937 mt(seed);
uniform_real_distribution<double> rand_dis(0.0, 1.0);

//-----
#define p 19937
#define q 7083//(x^p+x^q+1)
#define length 31


//y[] have randum numbers
int l=N*N+p;// size of y[i]
unsigned int y[22437];//N*N + p

double run=0;//
int set=0;//index of y[i]
//-------

ofstream outputfile("N=50.txt");
ofstream outputfile_2("N=50_bin.txt");

//for calulate specific_heat each bins
int BIN_TIME = TIME/bin;


int main() {

  double heat_cap=0;
  double temp = T_i;
  double beta;
  clock_t start,end;


  if(N%2==0){
    rand_set();
    set =p;
    start = clock();

    for(int i = 0; i <step; i++) {
      beta =1/temp;
      cout << temp << endl;
      work_even(temp, beta, heat_cap);
      outputfile << temp << " " << heat_cap << endl;
      temp += (T_f - T_i) / step;
    }
    end = clock();
    printf("It takes %.2f second\n",(double)(end-start)/CLOCKS_PER_SEC);
  }
  else{
    cout <<"odd N is unavailable"<<endl;
  }
  return 0;
}

void work_even(double & temp,double & beta,double &heat_cap){
  //energy
  double energy=0;
  double energy_now=0;
  double energy_sum=0;
  double ene_2=0;
  double ene_2_ave=0;

  //thermodynamics proparties for bin
  double heat_cap_bin=0;
  double energy_bin=0;
  double energy_sum_bin=0;
  double ene_2_bin=0;
  double ene_2_ave_bin=0;

  set = p;

  double r_a;
  double r_b;

  //energy_change_1
  r_a = exp(-8*J*beta);
  //energy_change_2
  r_b = exp(-4*J*beta);

  initial();

  //burn_in
  for (int i = 0; i < burn_in; i++) {
  monte_iti_even_non_energy(r_a,r_b);
  rand_reset();
  set = p;
}

  //initial state of energy
  energy_cal(energy_now);

//MCMC
for (int i=0; i < TIME; i++) {
  monte_iti_even(r_a,r_b,energy_now);

  energy_sum += energy_now;
  energy_sum_bin += energy_now;
  ene_2 += energy_now * energy_now;
  ene_2_bin += energy_now * energy_now;
  rand_reset();
  set = p;

  //calulation in bins
  if((i+1) % BIN_TIME == 0){
  energy_bin=1/(double)BIN_TIME*energy_sum_bin;
  ene_2_ave_bin=1/(double)BIN_TIME*ene_2_bin;

  //specific_heat
  heat_cap_bin=beta*beta * ( ene_2_ave_bin - energy_bin*energy_bin )/(N*N);

  outputfile_2 <<temp<<" "<<heat_cap_bin<<endl;
  energy_sum_bin=0;
  ene_2_bin=0;
}


}


energy=1/(double)TIME*energy_sum;
ene_2_ave=1/(double)TIME*ene_2;

heat_cap=beta*beta * ( ene_2_ave - energy*energy )/(N*N);

cout << "energy is "<<energy<<endl;
cout << "Heat capacity is　"<<heat_cap<< endl;
}

void initial(){
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      double b = rand_dis(mt);
      if (b > 0.5) {
        a[i][j] = 1;
      } else {
        a[i][j] = -1;
      }
    }
  }
}

void monte_iti_even_non_energy(double & r_a, double & r_b){
  int around = 0;
  double r_1 = 0;
  double r_2 = 0;

  double b=0;
  double N_half = N/2;


  for (int i = 0; i < N; i++) {
    int row_1 = (i - 1 + N) % N;
    int row_2 = (i + 1) % N;

    if(i % 2 == 0){
      for (int j = 0; j < N_half; j++) {
        around = a[row_1][2 * j] + a[row_2][2 * j] + a[i][(2 * j - 1 + N) % N] +  a[i][(2 * j + 1) % N];


        if(a[i][2 * j] > 0){//注目しているスピンが上向きの場合
          //rたちは、エネルギーが上がる時のr_a,r_b以外は必要ない
          if(around == 4){//周囲が全て上向き
            r_1 = r_a;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_1 > b) {
              a[i][2 * j] = -a[i][2 * j];//　　エネルギーが下がる方へは必ず変化。エネルギーが上がる方向へは確率的に変化。
            }
          }
          else if(around == 2){//周囲の3つが上向き
            r_1 = r_b;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_1 > b) {
              a[i][2 * j] = -a[i][2 * j];//　　エネルギーが下がる方へは必ず変化。エネルギーが上がる方向へは確率的に変化。
            }
          }
          else{
            //本来乱数は必要ないが、漸化式を進めるために発生させておく。この過程が必要ないような改良が望まれる。
            y[set] = y[set-p] ^ y[set-q];
            set +=1;

            a[i][2 * j] = -a[i][2 * j];
          }
        }

        else{//注目しているスピンが下向きの場合
          if(around == -2){//周囲の1つが上向き
            r_1 = r_b;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_1 > b) {
              a[i][2 * j] = -a[i][2 * j];//　　エネルギーが下がる方へは必ず変化。エネルギーが上がる方向へは確率的に変化。
            }
          }
          else if(around == -4){//周囲が全て下向き
            r_1 = r_a;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_1 > b) {
              a[i][2 * j] = -a[i][2 * j];//　　エネルギーが下がる方へは必ず変化。エネルギーが上がる方向へは確率的に変化。
            }
          }
          else{
            y[set] = y[set-p] ^ y[set-q];
            set +=1;

            a[i][2 * j] = -a[i][2 * j];
          }
        }
      }
    }
    else{
      for (int j = 0; j < N_half; j++) {

        //ここでj=1,3,5...のサイトを指定。
        around = a[row_1][2*j+1] + a[row_2][2*j+1] + a[i][(2*j + N) % N] + a[i][(2*j + 2) % N];


        if( a[i][2*j+1] > 0){//注目しているスピンが上向きの場合
          if(around == 4){//周囲が全て上向き
            r_2 = r_a;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_2 > b) {
              a[i][2*j+1] = -a[i][2*j+1];
            }
          }
          else if(around == 2){//周囲の3つが上向き
            r_2 = r_b;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_2 > b) {
              a[i][2*j+1] = -a[i][2*j+1];
            }
          }
          else{
            y[set] = y[set-p] ^ y[set-q];
            set +=1;

            a[i][2*j+1] = -a[i][2*j+1];
          }
        }

        else{//注目しているスピンが下向きの場合
          if(around == -2){//周囲の1つが上向き
            r_2 = r_b;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_2 > b) {
              a[i][2*j+1] = -a[i][2*j+1];
            }
          }
          else if(around == -4){//周囲が全て下向き
            r_2 = r_a;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_2 > b) {
              a[i][2*j+1] = -a[i][2*j+1];
            }
          }
          else{
            y[set] = y[set-p] ^ y[set-q];
            set +=1;

            a[i][2*j+1] = -a[i][2*j+1];
          }
        }
      }

    }
  }

  for (int i = 0; i < N; i++) {
    int row_1 = (i - 1 + N) % N;
    int row_2 = (i + 1) % N;
    if(i % 2 == 0){
      for (int j = 0; j < 0.5*N; j++) {

        //ここでj=1,3,5...のサイトを指定。
        around = a[row_1][2*j+1] + a[row_2][2*j+1] + a[i][(2*j + N) % N] + a[i][(2*j + 2) % N];

        if( a[i][2*j+1] > 0){//注目しているスピンが上向きの場合
          if(around == 4){//周囲が全て上向き
            r_2 = r_a;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_2 > b) {
              a[i][2*j+1] = -a[i][2*j+1];
            }
          }
          else if(around == 2){//周囲の3つが上向き
            r_2 = r_b;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_2 > b) {
              a[i][2*j+1] = -a[i][2*j+1];
            }
          }
          else{
            y[set] = y[set-p] ^ y[set-q];
            set +=1;
            a[i][2*j+1] = -a[i][2*j+1];
          }
        }

        else{//注目しているスピンが下向きの場合
          if(around == -2){//周囲の1つが上向き
            r_2 = r_b;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_2 > b) {
              a[i][2*j+1] = -a[i][2*j+1];
            }

          }
          else if(around == -4){//周囲が全て下向き
            r_2 = r_a;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_2 > b) {
              a[i][2*j+1] = -a[i][2*j+1];
            }
          }
          else{
            y[set] = y[set-p] ^ y[set-q];
            set +=1;

            a[i][2*j+1] = -a[i][2*j+1];
          }
        }
      }
    }
    else{
      for (int j = 0; j < N_half; j++) {

        //ここでj=0,2,4,...のサイトを指定。
        around = a[row_1][2 * j] + a[row_2][2 * j] + a[i][(2 * j - 1 + N) % N] + a[i][(2 * j + 1) % N];

        if(a[i][2 * j] > 0){//注目しているスピンが上向きの場合
          if(around == 4){//周囲が全て上向き
            r_1 = r_a;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_1 > b) {
              a[i][2 * j] = -a[i][2 * j];
            }
          }
          else if(around == 2){//周囲の3つが上向き
            r_1 = r_b;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_1 > b) {
              a[i][2 * j] = -a[i][2 * j];
            }
          }
          else{
            y[set] = y[set-p] ^ y[set-q];
            set +=1;
            a[i][2 * j] = -a[i][2 * j];
          }
        }
        else{//注目しているスピンが下向きの場合
          if(around == -2){//周囲の1つが上向き
            r_1 = r_b;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_1 > b) {
              a[i][2 * j] = -a[i][2 * j];
            }
          }
          else if(around == -4){//周囲が全て下向き
            r_1 = r_a;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_1 > b) {
              a[i][2 * j] = -a[i][2 * j];
            }
          }
          else{
            y[set] = y[set-p] ^ y[set-q];
            set +=1;
            a[i][2 * j] = -a[i][2 * j];
          }
        }
      }

    }
  }


}

void monte_iti_even(double & r_a, double & r_b,double & energy_now){
  //around は注目しているスピンの周囲のスピンの合計
  int around = 0;
  double r_1 = 0;
  double r_2 = 0;
  //１つフリップさせた時にエネルギーが上がるかどうかの指標がr。

  double b=0;//フリップに使う乱数
  double N_half = N/2;


  for (int i = 0; i < N; i++) {
    //配列の列の部分を先に計算しておく。
    int row_1 = (i - 1 + N) % N;
    int row_2 = (i + 1) % N;

    if(i % 2 == 0){
      for (int j = 0; j < N_half; j++) {
        around = a[row_1][2 * j] + a[row_2][2 * j] + a[i][(2 * j - 1 + N) % N] +  a[i][(2 * j + 1) % N];


        if(a[i][2 * j] > 0){//注目しているスピンが上向きの場合
          //rたちは、エネルギーが上がる時のr_a,r_b以外は必要ない
          if(around == 4){//周囲が全て上向き
            r_1 = r_a;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_1 > b) {
              a[i][2 * j] = -a[i][2 * j];//　　エネルギーが下がる方へは必ず変化。エネルギーが上がる方向へは確率的に変化。
              //エネルギーを+８J
              energy_now += 8*J;
            }
          }
          else if(around == 2){//周囲の3つが上向き
            r_1 = r_b;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_1 > b) {
              a[i][2 * j] = -a[i][2 * j];//　　エネルギーが下がる方へは必ず変化。エネルギーが上がる方向へは確率的に変化。
              //エネルギーを+4J
              energy_now += 4*J;
            }
          }
          else{
            //本来乱数は必要ないが、漸化式を進めるために発生させておく。この過程が必要ないような改良が望まれる。
            y[set] = y[set-p] ^ y[set-q];
            set +=1;

            a[i][2 * j] = -a[i][2 * j];
            //エネルギーは、J*around*2だけ変化する。around==0なら変化なし、around==-2なら-4J,around==-4なら-8Jされる。
            energy_now += 2*J*around;
          }
        }

        else{//注目しているスピンが下向きの場合
          if(around == -2){//周囲の1つが上向き
            r_1 = r_b;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_1 > b) {
              a[i][2 * j] = -a[i][2 * j];//　　エネルギーが下がる方へは必ず変化。エネルギーが上がる方向へは確率的に変化。
              //エネルギーを+4J
              energy_now += 4*J;
            }
          }
          else if(around == -4){//周囲が全て下向き
            r_1 = r_a;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_1 > b) {
              a[i][2 * j] = -a[i][2 * j];//　　エネルギーが下がる方へは必ず変化。エネルギーが上がる方向へは確率的に変化。
              //エネルギーを+８J
              energy_now += 8*J;
            }
          }
          else{
            y[set] = y[set-p] ^ y[set-q];
            set +=1;

            a[i][2 * j] = -a[i][2 * j];
            //エネルギーは、J*(-around)*2だけ変化する。around==0なら変化なし、around==2なら-4J,around==4なら-8Jされる。
            energy_now += -2*J*around;
          }
        }
      }
    }
    else{
      for (int j = 0; j < N_half; j++) {

        //ここでj=1,3,5...のサイトを指定。
        around = a[row_1][2*j+1] + a[row_2][2*j+1] + a[i][(2*j + N) % N] + a[i][(2*j + 2) % N];


        if( a[i][2*j+1] > 0){//注目しているスピンが上向きの場合
          if(around == 4){//周囲が全て上向き
            r_2 = r_a;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_2 > b) {
              a[i][2*j+1] = -a[i][2*j+1];
              energy_now += 8*J;
            }
          }
          else if(around == 2){//周囲の3つが上向き
            r_2 = r_b;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_2 > b) {
              a[i][2*j+1] = -a[i][2*j+1];
              energy_now += 4*J;
            }
          }
          else{
            y[set] = y[set-p] ^ y[set-q];
            set +=1;

            a[i][2*j+1] = -a[i][2*j+1];
            //エネルギーは、J*around*2だけ変化する。around==0なら変化なし、around==-2なら-4J,around==-4なら-8Jされる。
            energy_now += 2*J*around;
          }
        }

        else{//注目しているスピンが下向きの場合
          if(around == -2){//周囲の1つが上向き
            r_2 = r_b;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_2 > b) {
              a[i][2*j+1] = -a[i][2*j+1];
              energy_now += 4*J;
            }
          }
          else if(around == -4){//周囲が全て下向き
            r_2 = r_a;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_2 > b) {
              a[i][2*j+1] = -a[i][2*j+1];
              energy_now += 8*J;
            }
          }
          else{
            y[set] = y[set-p] ^ y[set-q];
            set +=1;

            a[i][2*j+1] = -a[i][2*j+1];
            //エネルギーは、J*(-around)*2だけ変化する。around==0なら変化なし、around==2なら-4J,around==4なら-8Jされる。
            energy_now += -2*J*around;
          }
        }
      }

    }
  }

  for (int i = 0; i < N; i++) {
    int row_1 = (i - 1 + N) % N;
    int row_2 = (i + 1) % N;
    if(i % 2 == 0){
      for (int j = 0; j < 0.5*N; j++) {

        //ここでj=1,3,5...のサイトを指定。
        around = a[row_1][2*j+1] + a[row_2][2*j+1] + a[i][(2*j + N) % N] + a[i][(2*j + 2) % N];

        if( a[i][2*j+1] > 0){//注目しているスピンが上向きの場合
          if(around == 4){//周囲が全て上向き
            r_2 = r_a;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_2 > b) {
              a[i][2*j+1] = -a[i][2*j+1];
              energy_now += 8*J;
            }
          }
          else if(around == 2){//周囲の3つが上向き
            r_2 = r_b;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_2 > b) {
              a[i][2*j+1] = -a[i][2*j+1];
              energy_now += 4*J;
            }
          }
          else{
            y[set] = y[set-p] ^ y[set-q];
            set +=1;
            a[i][2*j+1] = -a[i][2*j+1];
            //エネルギーは、J*around*2だけ変化する。around==0なら変化なし、around==-2なら-4J,around==-4なら-8Jされる。
            energy_now += 2*J*around;
          }
        }

        else{//注目しているスピンが下向きの場合
          if(around == -2){//周囲の1つが上向き
            r_2 = r_b;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_2 > b) {
              a[i][2*j+1] = -a[i][2*j+1];
              energy_now += 4*J;
            }

          }
          else if(around == -4){//周囲が全て下向き
            r_2 = r_a;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_2 > b) {
              a[i][2*j+1] = -a[i][2*j+1];
              energy_now += 8*J;
            }
          }
          else{
            y[set] = y[set-p] ^ y[set-q];
            set +=1;

            a[i][2*j+1] = -a[i][2*j+1];
            //エネルギーは、J*(-around)*2だけ変化する。around==0なら変化なし、around==2なら-4J,around==4なら-8Jされる。
            energy_now += -2*J*around;
          }
        }
      }
    }
    else{
      for (int j = 0; j < N_half; j++) {

        //ここでj=0,2,4,...のサイトを指定。
        around = a[row_1][2 * j] + a[row_2][2 * j] + a[i][(2 * j - 1 + N) % N] + a[i][(2 * j + 1) % N];

        if(a[i][2 * j] > 0){//注目しているスピンが上向きの場合
          if(around == 4){//周囲が全て上向き
            r_1 = r_a;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_1 > b) {
              a[i][2 * j] = -a[i][2 * j];
              energy_now += 8*J;
            }
          }
          else if(around == 2){//周囲の3つが上向き
            r_1 = r_b;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_1 > b) {
              a[i][2 * j] = -a[i][2 * j];
              energy_now += 4*J;
            }
          }
          else{
            y[set] = y[set-p] ^ y[set-q];
            set +=1;

            a[i][2 * j] = -a[i][2 * j];
            //エネルギーは、J*around*2だけ変化する。around==0なら変化なし、around==-2なら-4J,around==-4なら-8Jされる。
            energy_now += 2*J*around;
          }
        }
        else{//注目しているスピンが下向きの場合
          if(around == -2){//周囲の1つが上向き
            r_1 = r_b;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_1 > b) {
              a[i][2 * j] = -a[i][2 * j];
              energy_now += 4*J;
            }
          }
          else if(around == -4){//周囲が全て下向き
            r_1 = r_a;

            y[set] = y[set-p] ^ y[set-q];
            b = (double)y[set]/RAND_MAX;
            set +=1;

            if (r_1 > b) {
              a[i][2 * j] = -a[i][2 * j];
              energy_now += 8*J;
            }
          }
          else{
            y[set] = y[set-p] ^ y[set-q];
            set +=1;

            a[i][2 * j] = -a[i][2 * j];
            //エネルギーは、J*(-around)*2だけ変化する。around==0なら変化なし、around==2なら-4J,around==4なら-8Jされる。
            energy_now += -2*J*around;
          }
        }
      }

    }
  }
}

void energy_cal(double & energy_now){
  //一つのサイトが右と下へ二本づつ腕を伸ばしていると考えればOK
  for (int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      energy_now += -J * a[i][j]*(a[i][(j+1) % N ] + a[(i+1) % N][j]);
    }
  }
}


void rand_set(){
  //x[]は初期値の収容
  unsigned int x[p];
  unsigned int k=0;

  //y[i]の初期化
  for(int i=0;i<l;i++){
    y[i]=0;
  }

  //初期値の設定
  for(int i=0; i<p; i++){
    x[i]=rand();//0 ~ 2^31-1の間の整数をランダムでだす。
  }
  //cout <<"rand_set" <<endl;

  //漸化式の構築。
  //初期条件の代入
  for(int i=0; i<p;i++){
    y[i]=x[i];
  }
}

void rand_reset(){
  /*乱数をとり直す。
  初期条件はy[0]~y[p-1]であり、
  y[p]~y[p+N^2-1]が乱数として使われていた。
  次から乱数は
  y[p+N^2]=y[p+N^2-p]+y[p+N^2-q]
  から始まるため、
  y[N^2]~y[p+N^2-1]
  を新しい初期条件として代入して、それ以降は0に初期化する。

  setの取り直しもやる
  */


  //初期条件の取り直し
  for(int i=0;i<p;i++){
    y[i] = y[i + N*N];
  }

  for(int i=p;i<l;i++){
    y[i] = 0;
  }

}
