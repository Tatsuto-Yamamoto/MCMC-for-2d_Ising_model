//計算時間の計測。乱数のリセットとエネルギーの計算を削除している。
//bin数
//IF文の中を軽くなった。
//ising.cxxと同じ計算結果を与えることを確認済み。

# include <iostream>
# include <math.h>
# include <iomanip>
# include <random>
# include <fstream>
# include <time.h>
# include <stdio.h>
# include <limits.h>


using namespace std;
//std::coutのstd書かなくていい

//N,TIME,filename,l,T_iは手で直すこと。

#define N 50 //格子の数
#define TIME 100000000
#define step 6
#define seed 3
#define bin 10
#define T_i 2.25
#define T_f 2.32
#define burn_in 2000

#define J 1.0//結合定数(H=-Jx_ix_i+1)

int a[N][N];

//初期条件の設定
void initial();
//市松模様のモンテカルロ
void monte_iti_even(double &,double &,double &);
//市松模様のモンテカルロを、エネルギーの計算なしで行う(burn_inで用いる)
void monte_iti_even_non_energy(double &,double &);
//エネルギーの計算
void energy_cal(double &);
//全体の挙動
void work_even(double &,double &,double &);
//乱数の初期値を与える
void rand_set();
//乱数の初期化
void rand_reset();


//random_device rnd;
mt19937 mt(seed);
uniform_real_distribution<double> rand_dis(0.0, 1.0);

//-----
//乱数生成のための
#define p 19937
#define q 7083//(x^p+x^q+1)
#define length 31


//y[]は実際に作る乱数列
int l=N*N+p;// lはy[i]の配列のサイズ
unsigned int y[22437];//N*N + p

double run=0;//乱数
int set=0;//y[i]の引数。
//-------


//ファイルに出力するために使う関数。全体での比熱をファイルに出力する。
ofstream outputfile("N=50.txt");
//ファイルに出力するために使う関数。binごとの比熱をファイルに出力する。
ofstream outputfile_2("N=50_bin.txt");

//binごとに比熱を計算したいため定義する、あとで使う指数。workの外でも構わないが、一応ここに。
int BIN_TIME = TIME/bin;


int main() {

  double heat_cap=0;
  double temp = T_i;//温度
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
    printf("%.2f秒かかりました\n",(double)(end-start)/CLOCKS_PER_SEC);
  }
  else{
    cout <<"odd N is unavailable"<<endl;
  }
  return 0;
}

void work_even(double & temp,double & beta,double &heat_cap){
  //物理量の設定
  double energy=0;
  double energy_now=0;
  double energy_sum=0;
  double ene_2=0;
  double ene_2_ave=0;

  //binごとの物理量の計算
  double heat_cap_bin=0;
  double energy_bin=0;
  double energy_sum_bin=0;
  double ene_2_bin=0;
  double ene_2_ave_bin=0;


  set = p;//初期条件以外での乱数を

  //スピンのフリップを判定する際の指数の値
  double r_a;
  double r_b;


  //周囲の４つと、フリップ前の対象のスピンの向きが同じ場合,エネルギーは＋８J変化する。
  r_a = exp(-8*J*beta);
  //周囲の３つと、フリップ前の対象のスピンの向きが同じ場合,エネルギーは＋４J変化する。
  r_b = exp(-4*J*beta);



  //初期条件
  initial();

  //ここでburn_inを行う
  for (int i = 0; i < burn_in; i++) {
  monte_iti_even_non_energy(r_a,r_b);
  rand_reset();
  set = p;
}

  //エネルギーの初期状態を決める。
  energy_cal(energy_now);
//energy_nowが、モンテカルロステップ内でのエネルギーの二乗になる。


//市松模様のモンテカルロ。フリップが終わったところでenegy_nowにその配位でのエネルギーを代入してある。

//モンテカルロ(市松模様)
for (int i=0; i < TIME; i++) {
  //毎ステップごとにエネルギーを計算して、それを足しあげる。
  monte_iti_even(r_a,r_b,energy_now);
  //ここでエネルギーを代入。
  energy_sum += energy_now;
  energy_sum_bin += energy_now;
  ene_2 += energy_now * energy_now;
  ene_2_bin += energy_now * energy_now;
  rand_reset();
  set = p;

  //binごとに比熱を計算して出力する
  //そのためには、平均をとるために、BIN_TIME個のエネルギーを持ってきてそれの和を作って、適宜リセットする
  //energy_calの部分でbinごとでのエネルギーの和を用意する必要あり


  if((i+1) % BIN_TIME == 0){
  //BIN_TIME個のエネルギーで和を

  //エネルギーの平均値
  energy_bin=1/(double)BIN_TIME*energy_sum_bin;
  //エネルギーの二乗平均
  ene_2_ave_bin=1/(double)BIN_TIME*ene_2_bin;

  //比熱の計算, (二乗の期待値)−(期待値の二乗)
  heat_cap_bin=beta*beta * ( ene_2_ave_bin - energy_bin*energy_bin )/(N*N);
  //binごとの比熱は、全体とは別のファイルに出力しておく
  outputfile_2 <<temp<<" "<<heat_cap_bin<<endl;

  //binごとに別々で比熱を計算するため、ここでリセット
  energy_sum_bin=0;
  ene_2_bin=0;
}


}


//ここで物理量を出力
//エネルギーの平均値
energy=1/(double)TIME*energy_sum;
//エネルギーの二乗平均
ene_2_ave=1/(double)TIME*ene_2;

//比熱の計算, (二乗の期待値)−(期待値の二乗)
heat_cap=beta*beta * ( ene_2_ave - energy*energy )/(N*N);


cout << "energy is "<<energy<<endl;
cout << "Heat capacity is　"<<heat_cap<< endl;


}

void initial(){
  //初期状態をランダムに決定する。
  //  cout << "initial state is"<<endl;

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      double b = rand_dis(mt);
      if (b > 0.5) {
        a[i][j] = 1;
      } else {
        a[i][j] = -1;
      }
      //    cout <<setw(2)<<right<< a[i][j] << "  ";
    }
    //  cout << endl;
  }
  //cout << "and then" << endl;
}

void monte_iti_even_non_energy(double & r_a, double & r_b){
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
