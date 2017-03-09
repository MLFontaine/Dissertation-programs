/*  MF
 *  g++ lyap.cpp -lmpfr -lgmp -o lyap
 *  ./lyap
 */

#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>

// MACROS
#define DISP_DIG 5 //0 gives the default large amount of displayed digits with mpfr_out_str();
#define INIT_DELTA 0.0000000000000000000000000000000000000000000000000000000000000001 //ITERS depends on this roughly
#define ITERS 1000 //10000
#define TRANSIENT_SKIP 100 //100000
#define ATTRACTOR_POINTS 1000 //100000

#define PARAM_RESOLUTION 0.001


// Global variables
mpfr_t cusp, u1, u2, u3, u4, u5, u6, u7, v0, v1, v2, v3, v4, v5, midv, radv, lefv, rigv, temp, temp2, rand_Num, mydelta, mylog;
gmp_randstate_t r_state;

int first_itinerary1=0;
int first_itinerary2=0;

long int count_fails=0;

void the_func(mpfr_t& u, mpfr_t& v, int& itinerary);
void iterate_to_separate(mpfr_t& u, mpfr_t& v, mpfr_t& uu, mpfr_t& vv, int& itinerary1, int& itinerary2);

int main()
{
    // Flags and counters
    int itinerary1=0;
    int itinerary2=0;

    // Don't need to be global variables
    mpfr_t u, v, uu, vv, u_hold, v_hold, uu_hold, vv_hold, toavg;


    // Initialize global variables with 256 bits
    mpfr_inits2(256, cusp, u1, u2, u3, u4, u5, u6, u7, v0, v1, v2, v3, v4, v5, midv, radv, lefv, rigv, u, v, uu, vv, temp, temp2, mydelta, mylog, toavg, rand_Num, u_hold, v_hold, uu_hold, vv_hold, (mpfr_ptr) 0);


    // Fixed function parameters
    mpfr_set_d (cusp, 1.5, MPFR_RNDN); //make equal to u5 usually
    mpfr_set_d (u1, 0, MPFR_RNDN);
    mpfr_set_d (u2, 0, MPFR_RNDN);
    mpfr_set_d (u3, .5, MPFR_RNDN);
    mpfr_set_d (u4, 1, MPFR_RNDN);
    mpfr_set_d (u5, 1.5, MPFR_RNDN);
    mpfr_set_d (u6, .2, MPFR_RNDN); //u6 and u7 are the reinjected part, intially .2 and 1.2 gave period 4, .2 and .6 give complexity
    mpfr_set_d (u7, .6, MPFR_RNDN);
    mpfr_set_d (v0, 0, MPFR_RNDN);
    mpfr_set_d (v1, .333333, MPFR_RNDN);
    mpfr_set_d (v2, .666666666, MPFR_RNDN);
    mpfr_set_d (v3, 1, MPFR_RNDN);
    mpfr_set_d (v4, .4, MPFR_RNDN);
    mpfr_set_d (v5, .6, MPFR_RNDN);

    // Parameters to be calculated
    mpfr_set_d (midv, 0, MPFR_RNDN);
    mpfr_set_d (radv, 0, MPFR_RNDN);
    mpfr_set_d (lefv, 0, MPFR_RNDN);
    mpfr_set_d (rigv, 0, MPFR_RNDN);

    // midv = (v3+v0)/2;
    mpfr_add(midv, v3, v0, MPFR_RNDN);
    mpfr_div_d(midv, midv, 2.0, MPFR_RNDN);

    // u1 = -midv;
    mpfr_neg (u1, midv, MPFR_RNDN);

    // Set delta nbrhd with macro, initialize average holder
    mpfr_set_d (mydelta, INIT_DELTA, MPFR_RNDN);
    mpfr_set_d(toavg, 0, MPFR_RNDN);


    // Setup random stuff, 2 steps

    // 1. get an integer random seed using system utility
    // http://www.cs.yale.edu/homes/aspnes/pinewiki/C%282f%29Randomization.html

    unsigned int seed;
    FILE *f;
    f = fopen("/dev/urandom", "r"); //urandom less robust but quicker than random
    fread(&seed, sizeof(seed), 1, f);
    fclose(f);
    //printf("%u\n", seed);
    seed = 908435589; //use to fix the seed with some integer instead of using /dev/urandom/

    // 2. using the seed, generate random numbers held in r_state, randinit_default is Mersenne Twister
    // https://gmplib.org/list-archives/gmp-bugs/2008-March/000972.html

    gmp_randinit_default(r_state);
    gmp_randseed_ui(r_state, seed);

    // uncomment to generate and display a randum number example
    //mpfr_urandom(rand_Num, r_state, MPFR_RNDN);
    //mpfr_out_str(stdout, 10, DISP_DIG, rand_Num, MPFR_RNDN);


    // Initialize random values
    //////If u or v are within delta of the boundary, uu and vv could overflow... probably need to add a test

    //u
    mpfr_urandom(rand_Num, r_state, MPFR_RNDN);
    mpfr_mul_d(u, rand_Num, 1.5, MPFR_RNDN);

    //v
    mpfr_urandom(rand_Num, r_state, MPFR_RNDN);
    mpfr_set(v, rand_Num, MPFR_RNDN);

    //skip specified number of transients for u and v, select nearby uu and vv next when done
    for (int ii=0; ii<TRANSIENT_SKIP; ii++)
    {
        the_func(u, v, itinerary1);
    }


    //uu=u+(2*rand-1)*delta //spread [0,1] to [-1,1] in the parenthesis, so we have a square delta nbrhd
    mpfr_urandom(rand_Num, r_state, MPFR_RNDN);
    mpfr_mul_d(rand_Num, rand_Num, 2, MPFR_RNDN);
    mpfr_sub_d(rand_Num, rand_Num, 1, MPFR_RNDN);
    mpfr_mul(rand_Num, rand_Num, mydelta, MPFR_RNDN);
    mpfr_add(uu, u, rand_Num, MPFR_RNDN);

    //vv=v+(2*rand-1)*delta
    mpfr_urandom(rand_Num, r_state, MPFR_RNDN);
    mpfr_mul_d(rand_Num, rand_Num, 2, MPFR_RNDN);
    mpfr_sub_d(rand_Num, rand_Num, 1, MPFR_RNDN);
    mpfr_mul(rand_Num, rand_Num, mydelta, MPFR_RNDN);
    mpfr_add(vv, v, rand_Num, MPFR_RNDN);

    mpfr_set(u_hold, u, MPFR_RNDN);
    mpfr_set(v_hold, v, MPFR_RNDN);

    mpfr_set(uu_hold, uu, MPFR_RNDN);
    mpfr_set(vv_hold, vv, MPFR_RNDN);

    //----------------------------------
    // use this block to see if different paramter sets separate
    //----------------------------------

    FILE *ff;
    ff = fopen("param_exps.txt", "wb");

    double u6_init=0;
    double u7_init=0;

    long int num_params = 0;

    for(int i=1; u6_init<=1.4001; i++)
    {
        u7_init = u6_init+PARAM_RESOLUTION;

        for(int j=1; u7_init<=1.5001; j++)
        {
            mpfr_set_d (u6, u6_init, MPFR_RNDN);
            mpfr_set_d (u7, u7_init, MPFR_RNDN);

            mpfr_out_str (ff, 10, DISP_DIG, u6, MPFR_RNDN);
            fprintf(ff, ", ");
            mpfr_out_str (ff, 10, DISP_DIG, u7, MPFR_RNDN);
            fprintf(ff, ", ");

            ///
            //u
            mpfr_urandom(rand_Num, r_state, MPFR_RNDN);
            mpfr_mul_d(u, rand_Num, 1.5, MPFR_RNDN);

            //v
            mpfr_urandom(rand_Num, r_state, MPFR_RNDN);
            mpfr_set(v, rand_Num, MPFR_RNDN);

            //skip specified number of transients for u and v, select nearby uu and vv next when done
            for (int ii=0; ii<TRANSIENT_SKIP; ii++)
            {
                the_func(u, v, itinerary1);
            }

            //uu=u+(2*rand-1)*delta //spread [0,1] to [-1,1] in the parenthesis, so we have a square delta nbrhd
            mpfr_urandom(rand_Num, r_state, MPFR_RNDN);
            mpfr_mul_d(rand_Num, rand_Num, 2, MPFR_RNDN);
            mpfr_sub_d(rand_Num, rand_Num, 1, MPFR_RNDN);
            mpfr_mul(rand_Num, rand_Num, mydelta, MPFR_RNDN);
            mpfr_add(uu, u, rand_Num, MPFR_RNDN);

            //vv=v+(2*rand-1)*delta
            mpfr_urandom(rand_Num, r_state, MPFR_RNDN);
            mpfr_mul_d(rand_Num, rand_Num, 2, MPFR_RNDN);
            mpfr_sub_d(rand_Num, rand_Num, 1, MPFR_RNDN);
            mpfr_mul(rand_Num, rand_Num, mydelta, MPFR_RNDN);
            mpfr_add(vv, v, rand_Num, MPFR_RNDN);
            ///

            iterate_to_separate(u, v, uu, vv, itinerary1, itinerary2);

            mpfr_out_str(ff, 10, DISP_DIG, mylog, MPFR_RNDN);
            fprintf(ff, " \n");

            u7_init=u7_init+PARAM_RESOLUTION;
            num_params++;
        }
        u6_init=u6_init+PARAM_RESOLUTION;
    }

    printf("count_fails = %d", count_fails);
    putchar ('\n');
    printf("num_params = %d", num_params);
    putchar ('\n');
    printf("percent_fails = %f", (double)count_fails/num_params);
    putchar ('\n');

    fclose(ff);
    //----------------------------------
    // use this block to look at one fixed parameter set more closely
    //----------------------------------

    /*
    FILE *fp;
    fp = fopen("att_points.dat", "wb");


    for (int i=1; i<=ATTRACTOR_POINTS; i++)
    {
        mpfr_out_str (fp, 10, DISP_DIG, u, MPFR_RNDN);
        fprintf(fp, ", ");
        mpfr_out_str (fp, 10, DISP_DIG, v, MPFR_RNDN);
        fprintf(fp, " \n");

        iterate_to_separate(u, v, uu, vv, itinerary1, itinerary2);

        mpfr_set(u, u_hold, MPFR_RNDN);
        mpfr_set(v, v_hold, MPFR_RNDN);

        the_func(u, v, itinerary1);

        //should we skip this point if it starts on the left hand side?
        //if (itinerary1==3)
        //{
        //    continue;
        //}


        //uu=u+(2*rand-1)*delta //spread [0,1] to [-1,1] in the parenthesis, so we have a square delta nbrhd
        mpfr_urandom(rand_Num, r_state, MPFR_RNDN);
        mpfr_mul_d(rand_Num, rand_Num, 2, MPFR_RNDN);
        mpfr_sub_d(rand_Num, rand_Num, 1, MPFR_RNDN);
        mpfr_mul(rand_Num, rand_Num, mydelta, MPFR_RNDN);
        mpfr_add(uu, u, rand_Num, MPFR_RNDN);

        //vv=v+(2*rand-1)*delta
        mpfr_urandom(rand_Num, r_state, MPFR_RNDN);
        mpfr_mul_d(rand_Num, rand_Num, 2, MPFR_RNDN);
        mpfr_sub_d(rand_Num, rand_Num, 1, MPFR_RNDN);
        mpfr_mul(rand_Num, rand_Num, mydelta, MPFR_RNDN);
        mpfr_add(vv, v, rand_Num, MPFR_RNDN);

        mpfr_set(u_hold, u, MPFR_RNDN);
        mpfr_set(v_hold, v, MPFR_RNDN);

        mpfr_add(toavg, toavg, mylog, MPFR_RNDN);
    }

    putchar ('\n');
    mpfr_div_ui(toavg, toavg, ATTRACTOR_POINTS, MPFR_RNDN);
    printf("The average exponent is: ");
    mpfr_out_str (stdout, 10, DISP_DIG, toavg, MPFR_RNDN);
    printf("\n\n");

    fclose(fp);

    */

  return 0;
}


void iterate_to_separate(mpfr_t& u, mpfr_t& v, mpfr_t& uu, mpfr_t& vv, int& itinerary1, int& itinerary2)
{
    for (int i=1; i<=ITERS; i++)
    {
        // Iterate u,v and then uu,vv
        the_func(u, v, itinerary1);
        the_func(uu, vv, itinerary2);

        // Store the first itinerary position
        if (i==1)
        {
            first_itinerary1=itinerary1;
            first_itinerary2=itinerary2;

            /*
            printf("\n");
            printf("first_itinerary1= %d", first_itinerary1);
            putchar ('\n');
            printf("first_itinerary2= %d", first_itinerary2);
            putchar ('\n');
             */
        }

        // This means separation has been achieved
        if ( ((itinerary1==3 && itinerary2!=3) || (itinerary1!=3 && itinerary2==3)) ) //note that the itinerary tracks where it was before the current iteration
        {
            /*
            printf("u= ");
            mpfr_out_str (stdout, 10, DISP_DIG, u, MPFR_RNDN);
            putchar ('\n');
            printf("v= ");
            mpfr_out_str (stdout, 10, DISP_DIG, v, MPFR_RNDN);
            putchar ('\n');
            printf("itinerary1= %d", itinerary1);
            putchar ('\n');

            printf("uu= ");
            mpfr_out_str (stdout, 10, DISP_DIG, uu, MPFR_RNDN);
            putchar ('\n');
            printf("vv= ");
            mpfr_out_str (stdout, 10, DISP_DIG, vv, MPFR_RNDN);
            putchar ('\n');
            printf("itinerary2= %d", itinerary2);
            putchar ('\n');
            putchar ('\n');
             */

            //entropy calculation
            //mylog=-log(mydelta)/last_iteration;
            mpfr_log(mylog, mydelta, MPFR_RNDN);
            mpfr_div_ui(mylog, mylog, i, MPFR_RNDN);
            mpfr_neg(mylog, mylog, MPFR_RNDN);

            /*
            putchar ('\n');
            printf("mydelta= ");
            mpfr_out_str (stdout, 10, DISP_DIG, mydelta, MPFR_RNDN);
            putchar ('\n');

            printf("iterations= %d", i);
            putchar ('\n');

            printf("Exponent: ");
            mpfr_out_str(stdout, 10, DISP_DIG, mylog, MPFR_RNDN);
            printf("\n");
            //matlab effectively had (last_iter-1) originally, but I think that was wrong
             */


            break; //could also do i=ITERS if need to do anything after this conditional
        }
        else if (i==ITERS)
        {
            //printf("Failure to separate in %d iterations.\n", ITERS);
            count_fails++;
            mpfr_set_d(mylog, -1, MPFR_RNDN);
        }
    }
}



//make it so that the_func can be easily changed to different functions, probably wouldn't need itinerary return for others though

void the_func(mpfr_t& u, mpfr_t& v, int& itinerary)
{
    if (mpfr_lessequal_p(u1, u) && mpfr_less_p(u, u2))
    {
        //if u1<=u && u<u2 %LHS
        //x = (u6-u7)/(u1-u2).*(u-u1)+u6;
        //y = (v4-v5)/(v0-v3).*(v-v0)+v4;

        mpfr_sub(temp, u6, u7, MPFR_RNDN);
        mpfr_sub(temp2, u1, u2, MPFR_RNDN);
        mpfr_div(temp, temp, temp2, MPFR_RNDN);
        mpfr_sub(u, u, u1, MPFR_RNDN);
        mpfr_mul(u, temp, u, MPFR_RNDN);
        mpfr_add(u, u, u6, MPFR_RNDN);

        mpfr_sub(temp, v4, v5, MPFR_RNDN);
        mpfr_sub(temp2, v0, v3, MPFR_RNDN);
        mpfr_div(temp, temp, temp2, MPFR_RNDN);
        mpfr_sub(v, v, v0, MPFR_RNDN);
        mpfr_mul(v, temp, v, MPFR_RNDN);
        mpfr_add(v, v, v4, MPFR_RNDN);

        itinerary = 1;
    }
    else if(mpfr_lessequal_p(u2, u) && mpfr_lessequal_p(u, u3))
    {
        //u2<=u && u<=u3 %LRHS
        //x = (cusp-u2)/(u2-u3).*(u-u2)+cusp;
        //y = (v3-v2)/(v0-v3).*(v-v0)+v3;

        mpfr_sub(temp, cusp, u2, MPFR_RNDN);
        mpfr_sub(temp2, u2, u3, MPFR_RNDN);
        mpfr_div(temp, temp, temp2, MPFR_RNDN);
        mpfr_sub(u, u, u2, MPFR_RNDN);
        mpfr_mul(u, temp, u, MPFR_RNDN);
        mpfr_add(u, u, cusp, MPFR_RNDN);

        mpfr_sub(temp, v3, v2, MPFR_RNDN);
        mpfr_sub(temp2, v0, v3, MPFR_RNDN);
        mpfr_div(temp, temp, temp2, MPFR_RNDN);
        mpfr_sub(v, v, v0, MPFR_RNDN);
        mpfr_mul(v, temp, v, MPFR_RNDN);
        mpfr_add(v, v, v3, MPFR_RNDN);

        itinerary = 2;
    }
    else if(mpfr_less_p(u3, u) && mpfr_less_p(u, u4))
    {
        //u3<u && u<u4 %MRHS
        //y = (lefv-rigv)/(u3-u4).*(u-u3)+lefv; %used for x calculation
        //x = -( (radv).^2 - (y-midv).^2 ).^(1/2);

        //radv = midv-(v1-v0)/v3.*v;
        mpfr_sub(radv, v1, v0, MPFR_RNDN);
        mpfr_div(radv, radv, v3, MPFR_RNDN);
        mpfr_mul(radv, radv, v, MPFR_RNDN);
        mpfr_sub(radv, midv, radv, MPFR_RNDN);

        //lefv = midv + radv;
        mpfr_add(lefv, midv, radv, MPFR_RNDN);

        //rigv = midv - radv;
        mpfr_sub(rigv, midv, radv, MPFR_RNDN);

        //y
        mpfr_sub(temp, lefv, rigv, MPFR_RNDN);
        mpfr_sub(temp2, u3, u4, MPFR_RNDN);
        mpfr_div(temp, temp, temp2, MPFR_RNDN);
        mpfr_sub(u, u, u3, MPFR_RNDN);
        mpfr_mul(u, temp, u, MPFR_RNDN);
        mpfr_add(v, u, lefv, MPFR_RNDN);

        //x
        mpfr_sqr(temp, radv, MPFR_RNDN);
        mpfr_sub(temp2, v, midv, MPFR_RNDN);
        mpfr_sqr(temp2, temp2, MPFR_RNDN);
        mpfr_sub(u, temp, temp2, MPFR_RNDN);
        mpfr_sqrt(u, u, MPFR_RNDN);
        mpfr_neg(u, u, MPFR_RNDN);

        itinerary = 3;
    }
    else if(mpfr_lessequal_p(u4, u) && mpfr_lessequal_p(u, u5))
    {
        //u4<=u && u<=u5  %RRHS
        //x = (u2-u5)/(u4-u5).*(u-u4)+u2;
        //y = (v0-v1)/(v0-v3).*(v-v0);

        mpfr_sub(temp, u2, u5, MPFR_RNDN);
        mpfr_sub(temp2, u4, u5, MPFR_RNDN);
        mpfr_div(temp, temp, temp2, MPFR_RNDN);
        mpfr_sub(u, u, u4, MPFR_RNDN);
        mpfr_mul(u, temp, u, MPFR_RNDN);
        mpfr_add(u, u, u2, MPFR_RNDN);

        mpfr_sub(temp, v0, v1, MPFR_RNDN);
        mpfr_sub(temp2, v0, v3, MPFR_RNDN);
        mpfr_div(temp, temp, temp2, MPFR_RNDN);
        mpfr_sub(v, v, v0, MPFR_RNDN);
        mpfr_mul(v, temp, v, MPFR_RNDN);

        itinerary = 4;

    }
    else
        printf("OUT OF BOUNDS\n");
}






