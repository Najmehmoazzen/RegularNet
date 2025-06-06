
#include"Kuramoto.Version6.h"//import Internal library Kuramoto                                                            $$$$
#include <time.h>//import External library for calculate time                                                              $$$$
#include <iomanip>//                                                                                                       $$$$
int main(){                                                                     //@@@           Beginning main              ---

    const double* data=read_data("data.txt");                                   //@@@ read data from data.txt and write them---

    const int number_of_samples = int(data[8]);
    int sample=0 ;
    for (sample ;sample < number_of_samples ;sample++){ 
        cout << sample;
    //-------------------------------------------------------------------------------------------------------------------------
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    Read and definition data,          ---
    //@@@                                     data.txt and Example file          @@@@    Number_of_node,Phases_initial,     ---
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    frequency,adj,coupling,delay,time  ---

    const int Number_of_node = int(data[1]);                                    //@@@        N=Number_of_node=1000          ---
    cout << "|------------------------------------------------------|\n"<< endl;//@@@                                      ---
    const double* frequency_layer1 = read_1D_W("Layer1",Number_of_node);        //@@@        w=natural frequency      L1    ---
    //double* Phases_initial_layer1 = read_1D_I("origin1",Number_of_node);        //@@@        I=initial Phases         L1    ---
    double* Phases_initial_layer1 = generate_random_phases_and_save(Number_of_node,sample);

    const double* const* adj_layer1 = read_2D_A("Layer1",Number_of_node);       //@@@        A=adjacency matrix       L1    ---
    const double* const* Intrafrust_layer1 = read_2D_b("Layer1",Number_of_node);//@@@        b=Intralayer frustration L1    ---
    cout << "|------------------------------------------------------|\n"<< endl;//@@@                                      ---
    const int time_stationary = int(data[4] * 0.2);                             //@@@    example T=20 time_stationary= 10   ---
    const int Number_Steps_time_stationary = int(time_stationary / data[3]);    //@@@   for example T=20 dt=0.01 >> = 1000  ---
    int coupling_step = round(data[5]/data[6]);                                 //@@@                                       ---
    double* Phases_next_layer1 = new double[Number_of_node];                    //@@@    Definition Phases next             ---
    double* Phases_layer1_previous = for_loop_equal(Phases_initial_layer1);     //@@@               Phases changer          ---
    for (coupling_step;coupling_step <= int(data[7]/data[6]);coupling_step++){  //@@@                                       ---@
        double coupling=coupling_step*data[6];                                  //@@@            call coupling              ---@
        ostringstream ostrcoupling;                                             //@@@    declaring output string stream     ---@
        ostrcoupling << fixed << setprecision(2) << coupling;                   //@@@  Sending a number as a stream output  ---@
        string strcoupling = ostrcoupling.str();                                //@@@ the str() converts number into string ---@
        //ofstream Phases_layer1("Save/Phases(time)VS(Node)/L1_k="+               //@@@       create file for phases L1       ---@
        //                    strcoupling+"layer1.txt");                         //@@@                                       ---@
        double time_step = double(data[2]);                                     //@@@     reset time for new time           ---@  @        
        for (time_step;time_step < int(data[4]/data[3]);time_step++){           //@@@                                       ---@  @
            double time_loop=time_step*data[3];                                 //@@@                                       ---@  @
            Runge_Kutta_4(Number_of_node,                                       //@@@   Runge-Kutta 4th Order Method  L1    ---@  @
                        data[3],                                                //@@@                                       ---@  @
                        coupling,                                               //@@@                                       ---@  @
                        frequency_layer1,                                       //@@@                                       ---@  @
                        Intrafrust_layer1,                                      //@@@                                       ---@  @
                        adj_layer1,                                             //@@@                                       ---@  @
                        Phases_layer1_previous,                                 //@@@                                       ---@  @
                        Phases_next_layer1);                                    //@@@                                       ---@  @
            Phases_layer1_previous = for_loop_equal(Phases_next_layer1);        //@@@           Back to the future L1       ---@  @
            check_scale(Number_of_node,Phases_layer1_previous);                 //@@@       scale phases in -pi tp pi L1    ---@  @
            //Phases_layer1 << time_loop << '	';                                 //@@@                                       ---@  @
            //for (int i = 0; i < Number_of_node; i++){                           //@@@                                       ---@  @
            //    Phases_layer1 << std::fixed << std::setprecision(2) <<          //@@@                                       ---@  @
           //                                Phases_layer1_previous[i] << '	';   //@@@                                       ---@  @
           //}                                                                    //@@@                                       ---@  @
           //Phases_layer1 << endl;                                               //@@@                                       ---@  @
        }                                                                       //@@@                                       ---@  @
        cout<<"k=" <<strcoupling <<endl;                                        //@@@                                       ---@
        //Phases_layer1.close();                                                  //@@@                                       ---@
    }                                                                           //@@@                                       ---@
    //ofstream Last_Phase_layer1("Save/Last_Phase/layer1.txt");                   //@@@                                       ---
    //for (int i = 0; i < Number_of_node; i++){                                   //@@@                                       ---
    //    Last_Phase_layer1 << Phases_layer1_previous[i] << endl;                 //@@@--->   print last coupling phases      ---
    //}                                                                           //@@@                                       ---
    //Last_Phase_layer1.close();                                                  //@@@                                       ---

    ostringstream filename_stream;
    filename_stream << "Save/Last_Phase/Last_Phase_S_" << sample << ".txt";
    string filename = filename_stream.str();
    ofstream Last_Phase_layer1(filename);

    if (!Last_Phase_layer1) {
        std::cerr << "ERROR: Cannot open file " << filename << " for writing.\n";
    } 
    else {
        for (int i = 0; i < Number_of_node; i++) {
            Last_Phase_layer1 << Phases_layer1_previous[i] << endl;
        }
        Last_Phase_layer1.close();
        cout << "Last phase of sample " << sample << " saved to " << filename << "\n";
    }
}
    //delete Phases_layer1_previous;                                              //@@@                                       ---
    //delete Phases_next_layer1;                                                  //@@@                                       ---
        
    return 0;                                                                   //@@@     dont return any thing             ---
}                                                                               //@@@                                       ---
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@---
//-----------------------------------------------------------------------------------------------------------------------------
