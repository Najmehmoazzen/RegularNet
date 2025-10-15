
#include"Kuramoto.Version6.h"//import Internal library Kuramoto                                                            $$$$
#include <time.h>//import External library for calculate time                                                              $$$$
#include <iomanip>//                                                                                                       $$$$
int main(){                                                                     //@@@           Beginning main              ---
    const double* data=read_data("data.txt");                                   //@@@ read data from data.txt and write them---
    const int number_of_samples = int(data[8]);
    
    // %%%%%%%%%%%%%%%%%%%%%%%%% sample loop %%%%%%%%%%%%%%%%%%%%%
    int sample=0 ;
    for (sample ;sample < number_of_samples ;sample++){ 
        cout << sample;
        const int Number_of_node = int(data[1]);
        cout << "|------------------------------------------------------|\n"<< endl;
        double* Phases_initial_layer1 = generate_random_phases_and_save(Number_of_node,sample);
        const double* const* adj = read_2D_A("Layer1",Number_of_node);
        cout << "|------------------------------------------------------|\n"<< endl;
        const int time_stationary = int(data[4] * 0.2);
        const int Number_Steps_time_stationary = int(time_stationary / data[3]);
        double* Phases_next_layer1 = new double[Number_of_node];
        double* Phases_previous = for_loop_equal(Phases_initial_layer1);
        
        // %%%%%%%%%%%%%%%%%%%%%%%%% time loop %%%%%%%%%%%%%%%%%%%%%
        //ofstream Phases("Save/Phase/S"+to_string(sample)+".txt");
        //int counter_start=0;
        //int snapshot_phase_save=10
        ofstream synctime("Save/sync vs time/r_t_S" + to_string(sample) + ".txt");
        double time_step = double(data[2]);
        for (time_step;time_step < int(data[4]/data[3]);time_step++){
            double time_loop=time_step*data[3]; 
            Runge_Kutta_4(Number_of_node,data[3],adj,Phases_previous,Phases_next_layer1);  
            Phases_previous = for_loop_equal(Phases_next_layer1);
            check_scale(Number_of_node,Phases_previous);
            /*
            if (counter_start<=snapshot_phase_save){
            Phases << time_loop << '	';
            for (int i = 0; i < Number_of_node; i++){
                Phases << std::fixed << std::setprecision(2) <<Phases_previous[i] << '	';   
            }
            counter_start++;
            //Phases << endl;
            }
            */
            synctime<<order_parameter(Number_of_node,Phases_previous)<<endl;
        }
        synctime.close();
        //Phases.close();
        // _____________________ time loop _____________________
        save_lastphase(sample,Number_of_node,Phases_previous);
    }
    return 0;
}
