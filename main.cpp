// This is a 3-species May-Leonard Monte Carlo simulation using evolution
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <random>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
using namespace std;

ofstream out;
//Creating our random number generator
random_device rd{};
mt19937 rng{ rd() };
uniform_real_distribution<double> dist{ 0.0, 1.0 }; //This generates a random number between 0 and 1 Usage example: dist(rng)

double normalCDF(double value);
double NormalCDFInverse(double p);
double RationalApproximation(double t);

double RationalApproximation(double t)
{
    // Abramowitz and Stegun formula 26.2.23.
    // The absolute value of the error should be less than 4.5 e-4.
    double c[] = {2.515517, 0.802853, 0.010328};
    double d[] = {1.432788, 0.189269, 0.001308};
    return t - ((c[2]*t + c[1])*t + c[0]) /
               (((d[2]*t + d[1])*t + d[0])*t + 1.0);
}

double NormalCDFInverse(double p)
{
    if (p <= 0.0 || p >= 1.0)
    {
        std::stringstream os;
        os << "Invalid input argument (" << p
           << "); must be larger than 0 but less than 1.";
        throw std::invalid_argument( os.str() );
    }

    // See article above for explanation of this section.
    if (p < 0.5)
    {
        // F^-1(p) = - G^-1(p)
        return -RationalApproximation( sqrt(-2.0*log(p)) );
    }
    else
    {
        // F^-1(p) = G^-1(1-p)
        return RationalApproximation( sqrt(-2.0*log(1-p)) );
    }
}

int random_integer(int sides) { // returns values from 1 to sides
      int roll = dist(rng) * sides+1;
    if (roll==sides+1){
        roll=sides;
    }
    return roll;
}

double gaussdist(double mean, double deviation){
    double lower_bound=-2;
    double upper_bound=2;
    double Phi_lb = normalCDF(lower_bound);
    double Phi_ub = normalCDF(upper_bound);

    //generate uniform (0,1)
    double u = dist(rng);

    //inverse CDF
    double Phi_x = (Phi_ub-Phi_lb)*u+Phi_lb;
    double F_inv_u = NormalCDFInverse(Phi_x);
    double x = mean + deviation*F_inv_u; // change this
     //cout << x<<endl;
     return x;

}
double normalCDF(double value)
{
    return 0.5 * erfc(-value * M_SQRT1_2);
}

class Particle{
public:
    int seperate_deaths;
    double deviation = 0.01; // this is the standard deviation of the predation probability
    double lowest_deathrate = 0.90;
    // in order to trigger the constant population version, set all birth and death rates to 0
    double r_deathrate = 0.1;
    double r_birthrate = 0.1;
    double r_predrate = 0.35;
    double r_diffrate = 0.3;
    double r_swaprate = 0.15;
    double p_deathrate = 0.1;
    double p_birthrate = 0.1;
    double p_predrate = 0.55;
    double p_diffrate = 0.1;
    double p_swaprate = 0.15;
    double s_deathrate = 0.1;
    double s_birthrate = 0.1;
    double s_predrate = 0.70;
    double s_diffrate = 0.05;
    double s_swaprate = 0.05;
public:
    int tag;
    double deathprob; // DO NOT LET THE PARTICLES BECOME IMMORTAL
    double birthprob;
    double predprob;
    double diffprob;
    double swapprob;
    Particle()=default;
    explicit Particle(int x){
        if (x==1) { // rock
            deathprob = r_deathrate / (r_deathrate + r_birthrate + r_predrate + r_diffrate + r_swaprate);
            birthprob = r_birthrate / (r_deathrate + r_birthrate +r_predrate + r_diffrate + r_swaprate);
            predprob = r_predrate / (r_deathrate + r_birthrate + r_predrate + r_diffrate + r_swaprate);
            diffprob = r_diffrate / (r_deathrate + r_birthrate + r_predrate + r_diffrate + r_swaprate);
            swapprob = r_swaprate / (r_deathrate + r_birthrate + r_predrate + r_diffrate + r_swaprate);
            tag = 1;
        } else if (x==2) { //paper
            deathprob = p_deathrate / (p_deathrate + p_birthrate + p_predrate + p_diffrate + p_swaprate);
            birthprob = p_birthrate / (p_deathrate + p_birthrate + p_predrate + p_diffrate + p_swaprate);
            predprob = p_predrate / (p_deathrate + p_birthrate + p_predrate + p_diffrate + p_swaprate);
            diffprob = p_diffrate / (p_deathrate + p_birthrate + p_predrate + p_diffrate + p_swaprate);
            swapprob = p_swaprate / (p_deathrate + p_birthrate + p_predrate + p_diffrate + p_swaprate);
            tag = 2;
        } else if (x==3) { //scissors
            deathprob = s_deathrate / (s_deathrate + s_birthrate + s_predrate + s_diffrate + s_swaprate);
            birthprob = s_birthrate / (s_deathrate + s_birthrate + s_predrate + s_diffrate + s_swaprate);
            predprob = s_predrate / (s_deathrate + s_birthrate + s_predrate + s_diffrate + s_swaprate);
            diffprob = s_diffrate / (s_deathrate + s_birthrate + s_predrate + s_diffrate + s_swaprate);
            swapprob = s_swaprate / (s_deathrate + s_birthrate + s_predrate + s_diffrate + s_swaprate);
            tag = 3;
        }
        if (deathprob!=0){
            seperate_deaths = 1;
        } else{
            seperate_deaths=0;
        }
    }
    explicit Particle(double *parent_info){
        double parent_deathprob = parent_info[0];
        double parent_birthprob = parent_info[1];
        double parent_predprob = parent_info[2];
        double parent_diffprob = parent_info[3];
        double parent_swapprob = parent_info[4];
        double parent_tag = parent_info[5];
        if (parent_deathprob>lowest_deathrate){
            seperate_deaths = 1;
        } else{
            seperate_deaths=0;
        }
        tag = parent_tag;
        int pass = 1;
        while (pass==1) {
            predprob = gaussdist(parent_predprob, deviation);
            double difference = parent_predprob - predprob;
            if (seperate_deaths == 1) {
                deathprob = parent_deathprob + difference / 4;
                birthprob = parent_birthprob + difference / 4;
                diffprob = parent_diffprob + difference / 4;
                swapprob = parent_swapprob + difference / 4;
                if (predprob>1 || predprob < 0 || deathprob<lowest_deathrate || deathprob > 1){
                } else {
                    pass=0;
                }
            } else {
                diffprob = parent_diffprob + difference / 2;
                swapprob = parent_swapprob + difference / 2;
                deathprob=0;
                birthprob=0;
                if (predprob>1||predprob<0||diffprob<0||diffprob>1){
                } else{
                    pass=0;
                }
            }
        }
    }

    void React(Particle ***array, int *N, int size, int site_capacity, int *site){ // remember to always update N and the array so that the particles are on the lowest slices
        double action = dist(rng);
        if (action < deathprob) {
            //cout << "death"<<endl;
            Death(array, N, size, site_capacity, site);
        } else if (action < deathprob+birthprob){
            //cout << "birth"<<endl;
            Birth(array, N, size, site_capacity, site);
        } else if (action < deathprob+birthprob+predprob) {
            //cout << "predate"<<endl;
            Predate(array, N, size, site_capacity, site);
        } else if (action < deathprob+birthprob+predprob+diffprob){
            //cout << "diffuse"<<endl;
            Diffuse(array, N, size, site_capacity, site);
        } else {
            //cout << "swap"<<endl;
            Swap(array, N, size, site_capacity, site);
        }
    }
private:
    void Death(Particle ***array, int *N, int size, int site_capacity, int* site){
        int slice = site[0];
        int x = site[1];
        int y = site[2];
        int species = array[slice][x][y].tag; // gives a number 1-3
        N[species-1]=N[species-1]-1;
        array[slice][x][y] = Particle();
        int site2[2];
        site2[0]=x;
        site2[1]=y;
        ShuffleDown(array, site_capacity, site2);
    }

    void Birth(Particle ***array, int *N, int size, int site_capacity, int* site){
        int slice = site[0];
        int x = site[1];
        int y = site[2];
        int right = x == (size - 1) ? 0 : x + 1;
        /* ^This line checks if x is equal to size - 1. If it is, it means x is at the right boundary of the array,
         * so the value of right is set to 0. Otherwise, x is incremented by 1, and that value is assigned to right. */
        int left = x == 0 ? size - 1 : x - 1;
        int up = y == (size - 1) ? 0 : y + 1;
        int down = y == 0 ? size - 1 : y - 1;

        double parent_info[6];
        parent_info[0]=array[slice][x][y].deathprob;
        parent_info[1]=array[slice][x][y].birthprob;
        parent_info[2]=array[slice][x][y].predprob;
        parent_info[3]=array[slice][x][y].diffprob;
        parent_info[4]=array[slice][x][y].swapprob;
        parent_info[5]=array[slice][x][y].tag;
        int parent_tag = array[slice][x][y].tag;

        int blanks = 0;
        if (array[site_capacity-1][x][y].tag==0){
            blanks=blanks+1;
        }
        if (array[site_capacity-1][x][up].tag==0){
            blanks=blanks+1;
        }
        if (array[site_capacity-1][x][down].tag==0){
            blanks=blanks+1;
        }
        if (array[site_capacity-1][right][y].tag==0){
            blanks=blanks+1;
        }
        if (array[site_capacity-1][left][y].tag==0){
            blanks=blanks+1;
        }
        //cout << "blanks = " << blanks <<endl;
        if (blanks == 0 ){
        } else{
            int i =0;
            while (i==0) {
                //cout << "loop run"<<endl;
                int num = random_integer(5);
                switch (num) {
                    case 1:
                        if (array[site_capacity-1][x][y].tag==0){
                            //cout << "same site"<< endl;
                            i=1;
                            N[parent_tag-1]=N[parent_tag-1]+1;
                            array[site_capacity-1][x][y]=Particle(parent_info);
                            int site2[2];
                            site2[0]=x;
                            site2[1]=y;
                            ShuffleDown(array, site_capacity, site2);
                        }
                        break;
                    case 2:
                        if (array[site_capacity-1][right][y].tag==0){
                            //cout << "right";
                            i=1;
                            N[parent_tag-1]=N[parent_tag-1]+1;
                            array[site_capacity-1][right][y]=Particle(parent_info);
                            int site2[2];
                            site2[0]=right;
                            site2[1]=y;
                            ShuffleDown(array, site_capacity, site2);
                        }
                        break;
                    case 3:
                        if (array[site_capacity-1][left][y].tag==0){
                            i=1;
                            //cout << "left"<<endl;
                            N[parent_tag-1]=N[parent_tag-1]+1;
                            array[site_capacity-1][left][y]=Particle(parent_info);
                            int site2[2];
                            site2[0]=left;
                            site2[1]=y;
                            ShuffleDown(array, site_capacity, site2);
                        }
                        break;
                    case 4:
                        if (array[site_capacity-1][x][up].tag==0){
                            i=1;
                            //cout << "up"<<endl;
                            N[parent_tag-1]=N[parent_tag-1]+1;
                            array[site_capacity-1][x][up]=Particle(parent_info);
                            int site2[2];
                            site2[0]=x;
                            site2[1]=up;
                            ShuffleDown(array, site_capacity, site2);
                        }
                        break;
                    case 5:
                        if (array[site_capacity-1][x][down].tag==0){
                            i=1;
                            //cout << "down"<<endl;
                            N[parent_tag-1]=N[parent_tag-1]+1;
                            array[site_capacity-1][x][down]=Particle(parent_info);
                            int site2[2];
                            site2[0]=x;
                            site2[1]=down;
                            ShuffleDown(array, site_capacity, site2);
                        }
                        break;
                    default:
                        cout << "Default case reached in the birth code"<< endl;
                }
            }
        }
    }

    void Predate(Particle ***array, int *N, int size, int site_capacity, int* site){
        int slice = site[0];
        int x = site[1];
        int y = site[2];
        int right = x == (size - 1) ? 0 : x + 1;
        int left = x == 0 ? size - 1 : x - 1;
        int up = y == (size - 1) ? 0 : y + 1;
        int down = y == 0 ? size - 1 : y - 1;

        double parent_info[6];
        parent_info[0]=array[slice][x][y].deathprob;
        parent_info[1]=array[slice][x][y].birthprob;
        parent_info[2]=array[slice][x][y].predprob;
        parent_info[3]=array[slice][x][y].diffprob;
        parent_info[4]=array[slice][x][y].swapprob;
        parent_info[5]=array[slice][x][y].tag;
        int parent_tag = array[slice][x][y].tag;
        int prey_tag;
        if (parent_tag==1){
            prey_tag=3;
        } else {
            prey_tag=parent_tag-1;
        }

        int isprey = 0;
        for (int level=0;level<site_capacity;level++){
            if (array[level][x][y].tag==prey_tag){
                isprey=isprey+1;
            }
        }
        for (int level=0;level<site_capacity;level++){
            if (array[level][right][y].tag==prey_tag){
                isprey=isprey+1;
            }
        }
        for (int level=0;level<site_capacity;level++){
            if (array[level][left][y].tag==prey_tag){
                isprey=isprey+1;
            }
        }
        for (int level=0;level<site_capacity;level++){
            if (array[level][x][up].tag==prey_tag){
                isprey=isprey+1;
            }
        }
        for (int level=0;level<site_capacity;level++){
            if (array[level][x][down].tag==prey_tag){
                isprey=isprey+1;
            }
        }
        if (isprey!=0){
            int i = 0;
            while (i==0){
                int box = random_integer(5);
                switch (box) {
                    case 1:
                        for (int level=0;level<site_capacity;level++){
                            if (array[level][x][y].tag==prey_tag){
                                array[level][x][y]=Particle(parent_info);
                                N[parent_tag-1]=N[parent_tag-1]+1;
                                N[prey_tag-1]=N[prey_tag-1]-1;
                                level=site_capacity;
                                i=1;
                            }
                        }
                        break;
                    case 2:
                        for (int level=0;level<site_capacity;level++){
                            if (array[level][right][y].tag==prey_tag){
                                array[level][right][y]=Particle(parent_info);
                                N[parent_tag-1]=N[parent_tag-1]+1;
                                N[prey_tag-1]=N[prey_tag-1]-1;
                                level=site_capacity;
                                i=1;
                            }
                        }
                        break;
                    case 3:
                        for (int level=0;level<site_capacity;level++){
                            if (array[level][left][y].tag==prey_tag){
                                array[level][left][y]=Particle(parent_info);
                                N[parent_tag-1]=N[parent_tag-1]+1;
                                N[prey_tag-1]=N[prey_tag-1]-1;
                                level=site_capacity;
                                i=1;
                            }
                        }
                        break;
                    case 4:
                        for (int level=0;level<site_capacity;level++){
                            if (array[level][x][up].tag==prey_tag){
                                array[level][x][up]=Particle(parent_info);
                                N[parent_tag-1]=N[parent_tag-1]+1;
                                N[prey_tag-1]=N[prey_tag-1]-1;
                                level=site_capacity;
                                i=1;
                            }
                        }
                        break;
                    case 5:
                        for (int level=0;level<site_capacity;level++){
                            if (array[level][x][down].tag==prey_tag){
                                array[level][x][down]=Particle(parent_info);
                                N[parent_tag-1]=N[parent_tag-1]+1;
                                N[prey_tag-1]=N[prey_tag-1]-1;
                                level=site_capacity;
                                i=1;
                            }
                        }
                        break;
                    default:
                        cout << "default case reached in predation"<< endl;
                }
            }

        }
    }

    void Diffuse(Particle ***array, int *N, int size, int site_capacity, int* site){
        int slice = site[0];
        int x = site[1];
        int y = site[2];
        int right = x == (size - 1) ? 0 : x + 1;
        /* ^This line checks if x is equal to size - 1. If it is, it means x is at the right boundary of the array,
         * so the value of right is set to 0. Otherwise, x is incremented by 1, and that value is assigned to right. */
        int left = x == 0 ? size - 1 : x - 1;
        int up = y == (size - 1) ? 0 : y + 1;
        int down = y == 0 ? size - 1 : y - 1;

        int blanks = 0;
        if (array[site_capacity-1][x][up].tag==0){
            blanks=blanks+1;
        }
        if (array[site_capacity-1][x][down].tag==0){
            blanks=blanks+1;
        }
        if (array[site_capacity-1][right][y].tag==0){
            blanks=blanks+1;
        }
        if (array[site_capacity-1][left][y].tag==0){
            blanks=blanks+1;
        }
        if (blanks==0){
        } else{
            int i = 0;
            while (i==0){
                int direction = random_integer(4);
                switch (direction) {
                    case 1:
                        //right
                        if (array[site_capacity-1][right][y].tag==0){
                            i=1;
                            array[site_capacity-1][right][y]=array[slice][x][y];
                            array[slice][x][y]=Particle();
                            int site2[2];
                            site2[0]=right;
                            site2[1]=y;
                            ShuffleDown(array, site_capacity, site2);
                        }
                        break;
                    case 2:
                        //left
                        if (array[site_capacity-1][left][y].tag==0){
                            i=1;
                            array[site_capacity-1][left][y]=array[slice][x][y];
                            array[slice][x][y]=Particle();
                            int site2[2];
                            site2[0]=left;
                            site2[1]=y;
                            ShuffleDown(array, site_capacity, site2);
                        }
                        break;
                    case 3:
                        //up
                        if (array[site_capacity-1][x][up].tag==0){
                            i=1;
                            array[site_capacity-1][x][up]=array[slice][x][y];
                            array[slice][x][y]=Particle();
                            int site2[2];
                            site2[0]=x;
                            site2[1]=up;
                            ShuffleDown(array, site_capacity, site2);
                        }
                        break;
                    case 4:
                        //down
                        if (array[site_capacity-1][x][down].tag==0){
                            i=1;
                            array[site_capacity-1][x][down]=array[slice][x][y];
                            array[slice][x][y]=Particle();
                            int site2[2];
                            site2[0]=x;
                            site2[1]=down;
                            ShuffleDown(array, site_capacity, site2);
                        }
                        break;
                    default:
                        i=1;
                        cout << "Default code reached in diffusion"<<endl;
                        break;
                }

            }
        }
        int site3[2];
        site3[0]=x;
        site3[1]=y;
        ShuffleDown(array, site_capacity, site3);
    }

    void Swap(Particle ***array, int *N, int size, int site_capacity, int* site){
        int slice = site[0];
        int x = site[1];
        int y = site[2];
        int right = x == (size - 1) ? 0 : x + 1;
        /* ^This line checks if x is equal to size - 1. If it is, it means x is at the right boundary of the array,
         * so the value of right is set to 0. Otherwise, x is incremented by 1, and that value is assigned to right. */
        int left = x == 0 ? size - 1 : x - 1;
        int up = y == (size - 1) ? 0 : y + 1;
        int down = y == 0 ? size - 1 : y - 1;

        int counter = 0;
        for (int level = 0;level<site_capacity;level++){
            if (array[level][right][y].tag==0){
                counter=counter+1;
            }
        }
        for (int level = 0;level<site_capacity;level++){
            if (array[level][left][y].tag==0){
                counter=counter+1;
            }
        }
        for (int level = 0;level<site_capacity;level++){
            if (array[level][x][up].tag==0){
                counter=counter+1;
            }
        }
        for (int level = 0;level<site_capacity;level++){
            if (array[level][x][down].tag==0){
                counter=counter+1;
            }
        }
        if (counter==0){
        } else {
            int swap_num = random_integer(counter);
            int encountered = 0;
            for (int level = 0;level<site_capacity;level++){
                if (array[level][right][y].tag!=0){
                    encountered=encountered+1;
                    if (encountered==swap_num){
                        //swap
                        Particle storage = array[slice][x][y];
                        array[slice][x][y]=array[level][right][y];
                        array[level][right][y]=storage;
                        int site2[2];
                        site2[0]=right;
                        site2[1]=y;
                        ShuffleDown(array, site_capacity, site2);
                    }
                }
            }
            for (int level = 0;level<site_capacity;level++){
                if (array[level][left][y].tag!=0){
                    encountered=encountered+1;
                    if (encountered==swap_num){
                        //swap
                        Particle storage = array[slice][x][y];
                        array[slice][x][y]=array[level][left][y];
                        array[level][left][y]=storage;
                        int site2[2];
                        site2[0]=left;
                        site2[1]=y;
                        ShuffleDown(array, site_capacity, site2);
                    }
                }
            }
            for (int level = 0;level<site_capacity;level++){
                if (array[level][x][up].tag!=0){
                    encountered=encountered+1;
                    if (encountered==swap_num){
                        //swap
                        Particle storage = array[slice][x][y];
                        array[slice][x][y]=array[level][x][up];
                        array[level][x][up]=storage;
                        int site2[2];
                        site2[0]=x;
                        site2[1]=up;
                        ShuffleDown(array, site_capacity, site2);
                    }
                }
            }
            for (int level = 0;level<site_capacity;level++){
                if (array[level][x][down].tag!=0){
                    encountered=encountered+1;
                    if (encountered==swap_num){
                        //swap
                        Particle storage = array[slice][x][y];
                        array[slice][x][y]=array[level][x][down];
                        array[level][x][down]=storage;
                        int site2[2];
                        site2[0]=x;
                        site2[1]=down;
                        ShuffleDown(array, site_capacity, site2);
                    }
                }
            }

        }
        int site3[2];
        site3[0]=x;
        site3[1]=y;
        ShuffleDown(array, site_capacity, site3);
    }

    void ShuffleDown(Particle ***array, int site_capacity, int* site){
        int x = site[0];
        int y = site[1];
        int list[site_capacity];
        for (int level = 0 ; level < site_capacity; level++){
            if (array[level][x][y].tag != 0){
                list[level]=1;
            } else {
                list[level] = 0;
            }
        }
        int num_placed=0;
        for (int level = 0; level < site_capacity; level++){
            if (list[level]==1){
                num_placed = num_placed + 1;
                if (num_placed-1!=level) {
                    array[num_placed - 1][x][y] = array[level][x][y];
                    array[level][x][y] = Particle();
                }
            }
        }
    }
};

void initialize_array(Particle ***array,int* N, int size,int site_capacity);
void update_array(Particle ***array, int *N, int size, int site_capacity);
int* pick_particle(Particle ***array, int *N, int size, int site_capacity);
double* get_average(Particle ***array, int *N, int size, int site_capacity);


int main() {
    out.open("DATA.txt");
    //Assumptions: square lattice, homogeneous
    const int size = 100;
    const int MC_steps = 2000;
    const int site_capacity = 3;
    const double r_density = 1;
    const double p_density = 1.2;
    const double s_density = 0.5;
    if (r_density+p_density+s_density > site_capacity){
      cout << "The total initial density is too high. " << endl;
      return 1;
    };
    int N[3];
    N[0] = r_density*size*size;
    N[1] = p_density*size*size;
    N[2] = s_density*size*size;
    Particle ***array = new Particle **[site_capacity];
    for (int slice=0; slice < site_capacity; slice++){
        array[slice] = new Particle * [size];
        for (int x = 0; x < size ; x++){
            array[slice][x] = new Particle[size];
            for (int y=0; y < size; y++){
                // initializes each particle object using the constructor
                array[slice][x][y] = Particle();
            }
        }
    }
    initialize_array(array, N, size, site_capacity);
    for (int t=0; t< MC_steps; t++){
       //out << (double) N[0] / (size*size) << " "<< (double) N[1] / (size*size) << " " << (double) N[2] / (size*size) << endl;
        double *pred_avg = get_average(array, N, size, site_capacity);
        //cout << "Pred_Avg's: "<<
        out << (double) pred_avg[0] << " " << (double) pred_avg[1] << " " << (double) pred_avg[2] <<" "<< (double) N[0] / (size*size) << " "<< (double) N[1] / (size*size) << " " << (double) N[2] / (size*size) <<endl;
        update_array(array, N, size, site_capacity);

    }
    out.close();
    return 0;
}

void update_array(Particle ***array, int *N, int size, int site_capacity) {
    int moves = N[0]+N[1]+N[2];
    for (int i=1; i <=moves; i++) {
        if (N[0]+N[1]+N[2]==0){
            return;
        }
        int* site = pick_particle(array, N, size, site_capacity);
        if (site==nullptr){
        } else {
            int slice = site[0];
            int x = site[1];
            int y = site[2];
            //cout << array[slice][x][y].deathprob << " " << array[slice][x][y].birthprob<<endl;
            //cout << array[slice][x][y].deathprob+array[slice][x][y].birthprob+array[slice][x][y].predprob+array[slice][x][y].diffprob+array[slice][x][y].swapprob<<endl;
            array[slice][x][y].React(array, N, size, site_capacity, site);
        }
    }
}

int* pick_particle(Particle ***array, int *N, int size, int site_capacity) {
    int *site = new int[3];
    int total = N[0]+N[1]+N[2];
    int roll = random_integer(total);
    int counter = 0;
    for (int x=0; x < size; x++){
        for (int y=0; y< size; y++){
            for (int slice = 0; slice < site_capacity; slice++) {
                if (array[slice][x][y].tag!=0){

                    counter=counter+1;
                    if (counter==roll){
                        site[0]=slice;
                        site[1]=x;
                        site[2]=y;
                        return site;
                    }
                }
            }
        }
    }
    cout << "A particle was not picked";
    return nullptr;
}


void initialize_array(Particle ***array, int*N, int size, int site_capacity){
    int num_r=N[0];
    int num_p=N[1];
    int num_s=N[2];
    for (int i = 0; i < size; i++){
        for (int j=0; j < size; j++){
            for (int k = 0; k<site_capacity; k++){
                array[k][i][j] = Particle();
            }
        }
    }
    for (int i = 1; i <=num_r; i++){ // places rocks
        int placed = 0;
        while (placed == 0) {
            int row = random_integer(size) - 1;
            int column = random_integer(size) - 1;
            for (int j = 0; j < site_capacity; j++) {
                if (placed==0) {
                    if (array[j][row][column].tag == 0) {
                        array[j][row][column]=Particle(1);
                        placed = 1;
                    }
                }
            }
        }
    }
    for (int i = 1; i <=num_p; i++){ // places paper
        int placed = 0;
        while (placed == 0) {
            int row = random_integer(size) - 1;
            int column = random_integer(size) - 1;
            for (int j = 0; j < site_capacity; j++) {
                if (placed==0) {
                    if (array[j][row][column].tag == 0) {
                        array[j][row][column]=Particle(2);
                        placed = 1;
                    }
                }
            }
        }
    }
    for (int i = 1; i <=num_s; i++){ // places scissors
        int placed = 0;
        while (placed == 0) {
            int row = random_integer(size) - 1;
            int column = random_integer(size) - 1;
            for (int j = 0; j < site_capacity; j++) {
                if (placed==0) {
                    if (array[j][row][column].tag == 0) {
                        array[j][row][column]=Particle(3);
                        placed = 1;
                    }
                }
            }
        }
    }
}

double* get_average(Particle ***array, int *N, int size, int site_capacity){
    static double averages[3];
    int num1=0;
    int num2 = 0;
    int num3=0;
    double pred1=0;
    double pred2=0;
    double pred3=0;
    for (int x=0; x<size;x++){
        for (int y=0; y < size; y++){
            for (int slice=0; slice<site_capacity;slice++){
                int tag = array[slice][x][y].tag;
                switch (tag) {
                    case 1: //rock
                        num1=num1+1;
                        pred1 = pred1 + array[slice][x][y].predprob;
                        break;
                    case 2:
                        num2=num2+1;
                        pred2 = pred2 + array[slice][x][y].predprob;
                        break;
                    case 3:
                        num3=num3+1;
                        pred3 = pred3 + array[slice][x][y].predprob;
                        break;
                    default:
                        break;
                }


            }
        }
    }
    if (num1==0){
        averages[0] = 0;
    } else {
        averages[0] = pred1 / num1;
    }
    if (num2==0){
        averages[1]=0;
    } else {
        averages[1] = pred2 / num2;
    }
    if (num3==0){
        averages[2]=0;
    } else {
        averages[2] = pred3/num3;
    }
    return averages;
}
