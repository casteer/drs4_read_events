/*
   Name:           read_binary.cpp
   Created by:     Stefan Ritt <stefan.ritt@psi.ch>
   Date:           July 30th, 2014

   Purpose:        Example file to read binary data saved by DRSOsc.
 
   Compile and run it with:
 
      gcc -o read_binary read_binary.cpp -lm -lstdc++
 
      ./read_binary <filename>

   This program assumes that a pulse from a signal generator is split
   and fed into channels #1 and #2. It then calculates the time difference
   between these two pulses to show the performance of the DRS board
   for time measurements.

   $Id: read_binary.cpp 21438 2014-07-30 15:00:17Z ritt $
*/

#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <cstdlib>
#include <math.h>
#include <vector>
#include <map>


#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/density.hpp>
#include <boost/accumulators/statistics/stats.hpp>
 
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>


typedef struct {
   char           time_header[4];
   char           bn[2];
   unsigned short board_serial_number;
} THEADER;

typedef struct {
   char           event_header[4];
   unsigned int   event_serial_number;
   unsigned short year;
   unsigned short month;
   unsigned short day;
   unsigned short hour;
   unsigned short minute;
   unsigned short second;
   unsigned short millisecond;
   unsigned short reserved1;
   char           bs[2];
   unsigned short board_serial_number;
   char           tc[2];
   unsigned short trigger_cell;
} EHEADER;

/*-----------------------------------------------------------------------------*/

class drs4_event
{

 public: 
	int ID;// Unique ID

	std::map<int, std::vector<double> > times;
	std::map<int, std::vector<double> > heights;
	std::map<int, int> hits;// number of hits in this event in each detector

	// Derived quantities
	std::map<int, std::vector<double> > inter_times;// time differences in same detector different pulses
	std::map<int, std::vector<double> > delay_times;// time differences to different detectors first pulse only 

	double lower_threshold;
	double upper_threshold;
	
	// from the header
	EHEADER event_header;


	// Constructor
	drs4_event(EHEADER eh, int new_id, double low_thresh, double upper_thresh){

		// Create a new event from the header 
		ID=new_id;
		lower_threshold=low_thresh;
		upper_threshold=upper_thresh;
		std::vector<double> v;
		for (int d=0 ; d<4 ; d++){
			hits.insert(std::pair<int,int>(d,0));		
			times.insert(std::pair<int,std::vector<double> >(d,v));		
			inter_times.insert(std::pair<int,std::vector<double> >(d,v));		
			heights.insert(std::pair<int,std::vector<double> >(d,v));		
		}


		// Save time from the event header
		event_header=eh;
	}





	double time_range;
		
	// function which fills event structure here
	void fill(double waveform[4][1024],double time[4][1024]){
  		
		double pulse_height=0;
		time_range=time[0][1022];


		for (int d=0 ; d<4 ; d++){
			double event_time;


			// may get half a pulse here?

			for (int i=0 ; i<1022 ; i++){

				// Waveform crosses threshold so we have a hit
         			if (waveform[d][i] <= lower_threshold && waveform[d][i+1] >= lower_threshold) {
					if(pulse_height<upper_threshold){
  	        				event_time = ((lower_threshold-waveform[d][i])/(waveform[d][i+1]-waveform[d][i]))*(time[d][i+1]-time[d][i])+time[d][i];
					}
         			}

		        	if (waveform[d][i] > pulse_height) pulse_height=waveform[d][i];

				// Finished processing the pulse
         			if (waveform[d][i] >= lower_threshold && waveform[d][i+1] <= lower_threshold) {
						// Save the event pulse height					
						hits[d]++;
	
						// Save the event time
						times[d].push_back(event_time);
						event_time=0;
						heights[d].push_back(pulse_height);
						pulse_height=0;	
				}
	 		}// samples loop

			// may get half a pulse here
			if (waveform[d][1022] >= lower_threshold) {
					// Save the event pulse height					
					hits[d]++;
	
					// Save the event time
					times[d].push_back(event_time);
					event_time=0;
					heights[d].push_back(pulse_height);
					pulse_height=0;	
			}



		}// channel loop
	}

	void fill_derived(){
			for (int d=0 ; d<4 ; d++){
				for (int h=0 ; h<hits[d]-1 ; h++){
					inter_times[d].push_back(times[d].at(h+1)-times[d].at(h));					
				}
			}
			if((hits[0]>0)&&(hits[1]>0)){			
				for (int h=0 ; h<hits[1] ; h++){
					delay_times[0].push_back(times[1].at(h)-times[0].at(0));						
				}
			}	
	}

	void print_times(){
		for (int d=0 ; d<4 ; d++){
			for (int h=0 ; h<hits[d] ; h++){
				printf("%f\n",times[d].at(h));
			}
		}
	}

	void print_times(int d){
			for (int h=0 ; h<hits[d] ; h++){
				printf("%f\n",times[d].at(h));
			}
	}

	void print_delay_times(){
			if((hits[0]>0)&&(hits[1]>0)){			
				for (int h=0 ; h<hits[1] ; h++){
					printf("%f\n",delay_times[0].at(h));
				}
			}
	}



	void print_inter_times(){
		for (int d=0 ; d<4 ; d++){
			for (int h=0 ; h<hits[d]-1 ; h++){
				printf("%f\n",inter_times[d].at(h));
			}
		}
	}
	void print_inter_times(int d){
			for (int h=0 ; h<hits[d]-1 ; h++){
				printf("%f\n",inter_times[d].at(h));
			}
	}

	void print_seq_times(){
		if((hits[0]>3)&&(hits[1]>3)) printf("%f\n",times[1].at(1)-times[0].at(0));
	}

	void print_multiplicity(){
		for (int d=0 ; d<4 ; d++){
			printf("%d\n",hits[d]);
		}
	}
	void print_multiplicity(int d){
			printf("%d\n",hits[d]);
	}


	void print_pulse_heights(){
		for (int d=0 ; d<4 ; d++){
			for (int h=0 ; h<hits[d]-1 ; h++){
				printf("%f\n",heights[d].at(h));
			}
		}
	}

	void print_pulse_heights(int d){
			for (int h=0 ; h<hits[d]-1 ; h++){
				printf("%f\n",heights[d].at(h));
			}
	}


	void print_event(){	
		for (int d=0 ; d<4 ; d++){
			std::cout << "# " << ID << " " << d << " " << hits[d] << std::endl;
			for (int h=0 ; h<hits[d] ; h++){
				std::cout << " " << times[d].at(h) << " " << heights[d].at(h) << std::endl;
			}
		}
	}



	void print(){
		printf("Event ID %d \n", ID);
		for (int d=0 ; d<4 ; d++){
			printf(" - Channel %d\n",d);
			printf(" - - Hits %d\n",hits[d]);
			for (int h=0 ; h<hits[d] ; h++){
				printf(" - - - Time %f\n",times[d].at(h));
				printf(" - - - Height %f\n",heights[d].at(h));
			}
		}
	}

};

/*
class drs4_cut
{

public: 

	double low_height;
	double high_height;

	int low_time;// time within the event
	int high_time;// time within the event

	int low_multiplicity;// time within the event
	int high_multiplicity;// time within the event

	drs4_cut(){
		low_height=-1000;
		high_height=1000;

		low_channel=0;
		high_channel=3;

		low_time=0;
		high_time=10000;
	}

	bool cut_event(event e){
		
		event e_new; 

		bool test_multiplicity = (e.hits[d]>low_multiplicity)&&(e.hits[d]<high_multiplicity);
		// Remove events outside the mult range 
		for (int d=0; d<4 ; d++){
			e.hits	



		bool test_time = (e.hits[d]>low_time)&&(e.hits[d]<high_time);


	}

}*/


class drs4_run
{

 public: 
	// Board that was used in the data-taking
	THEADER board_header;
	
	std::vector<drs4_event> events;

	// Return the average height? times

	std::vector<double> event_times;
	double time_range;
	
	drs4_run(THEADER th){
		board_header=th;
	};


	double add_event(drs4_event e){		
		events.push_back(e);
		
		event_times.push_back( (((double) e.event_header.hour)*60.0*60.0)+(((double) e.event_header.minute)*60.0)+e.event_header.second+((double) e.event_header.millisecond/1000.0));

		time_range=e.time_range;

	}

	double event_rate(){
		// Assumes insertion in chronological order
		return ((double) events.size())/(event_times.back()-event_times.front());
	}
	void print_event_times(){
		for (std::vector<double>::iterator it = event_times.begin() ; it != event_times.end(); ++it)
			std::cout << (*it) << std::endl;
	}

	void print_event_number(){
			std::cout << events.size() << std::endl;
	}	

	void print_event_rate(){
			std::cout << event_rate() << std::endl;
	}	

	void print_event_timescale(){
			std::cout << "# Time extent : " << time_range << std::endl;
	}	

};


using namespace boost;
using namespace boost::accumulators;
using namespace boost::program_options; 
 
typedef accumulator_set<double, features<tag::density,tag::count> > acc;
typedef iterator_range<std::vector<std::pair<double, double> >::iterator > histogram_type; 
 
template <typename T>
class data_filler 
{  
public:
  data_filler(){}
  T operator()() { return rand()/(T)RAND_MAX; }
};
 


int main(int argc, const char * argv[])
{
   THEADER th;
   EHEADER eh;
   char hdr[4];
   unsigned short voltage[1024];
   double waveform[4][1024], time[4][1024];
   float bin_width[4][1024];
   char rootfile[256];
   int i, j, ch, n, chn_index;
   double t1, t2, dt;
   char filename[256],h0_filename[256],h1_filename[256],t0_filename[256],t1_filename[256],t2_filename[256],hits0_filename[256],hits1_filename[256];

   int ndt;
   double lower_threshold, upper_threshold, sumdt, sumdt2;
   double pulse_height0=0.0;
   double pulse_height1=0.0;

	int flag=0;// whether above threshold or not
	int nhits1=0;
	int nhits2=0;

	options_description desc("DRS4 Data Analysis Program, the allowed options are");
	desc.add_options()
    	("help,h", "produce help message")
    	("channel", value<int>()->default_value(-1), "selects which channel to extract the data from")
    	("low_threshold", value<double>()->default_value(0.01), "event trigger threshold (V) ")
    	("high_threshold", value<double>()->default_value(1000.0), "removes pulses heights less than this threshold (V) ")
    	("pulse-times", value<bool>()->default_value(false), "output raw times data")
    	("pulse-seq_times", value<bool>()->default_value(false), "output sequenced times data = Time(hit 2, det1)-Time(hit 1,det0)")
    	("pulse-delay", value<bool>()->default_value(false), "output delay ch1-ch0 data")
    	("pulse-interevent", value<bool>()->default_value(false), "output interevent times")
    	("pulse-heights", value<bool>()->default_value(false),"output pulse height data")
    	("pulse-hits", value<bool>()->default_value(false),"output number of hits per event")
    	("run-times", value<bool>()->default_value(false),"output the event times")
    	("run-rate", value<bool>()->default_value(false),"output the event rate")
    	("run-number", value<bool>()->default_value(false),"output the number of events")
    	("start-event", value<int>()->default_value(0),"first event index to be processed")
    	("end-event", value<int>()->default_value(100000000),"last event to be processed")
    	("input,i", value<std::string>(),"input data file");

	variables_map vm;
	store(parse_command_line(argc, argv, desc), vm);
	notify(vm);    
	if (vm.count("help")) {
    		std::cout << desc << "\n";
    		return 1;
	}

	if (vm["pulse-times"].as<bool>()) {
		std::cout << "# Output : Raw Times" << std::endl;
	}

	if (vm["pulse-interevent"].as<bool>()) {
		std::cout << "# Output : Inter Event Times" << std::endl;
	}

	if (vm["pulse-heights"].as<bool>()) {
		std::cout << "# Output : Pulse heights " << std::endl;
	}

	lower_threshold=0.01;
	upper_threshold=1000.0;
	if (vm.count("threshold")) {
		lower_threshold=vm["low_threshold"].as<double>();
		upper_threshold=vm["high_threshold"].as<double>();
	}

	int channel=-1;
	if (vm.count("channel")) {
		channel=vm["channel"].as<int>();
//		std::cout << "# Output : Channel " << channel << std::endl;
	}

	if (!vm.count("input")) {
    		std::cout << "Error : Must inclue an input file" << std::endl;
    		std::cout << desc << "\n";
    		return 1;
	}else{

		std::string str_filename(vm["input"].as<std::string>());	
		strcpy(filename,str_filename.c_str());		
	}

	/*
		if (!vm.count("input")) {
		std::cin >> line;
		while(!line.empty()){
			data.push_back(atof(line.c_str()));
			line.clear();
			std::cin>>line;
			// std::cout << data.size() << std::endl;
		}
	}else{

		std::string filename(vm["input"].as<std::string>());	
		std::ifstream myfile(filename.c_str());// input filename
 		 if(myfile.is_open()) {
		 	while ( getline (myfile,line) ) {
 				data.push_back(atof(line.c_str()));
  			}
 		}
 		myfile.close();
	}
	*/



/*
   if (argc > 2){
	 threshold=atof(argv[2]);
   	printf("Threshold = %f\n",threshold);
      	strcpy(filename, argv[1]);
}  else {
      printf("Usage: read_binary <filename>\n");
      return 0;
   }



	// Remove
   	strcpy(t0_filename,"times0.dat");
   	strcpy(t1_filename,"times1.dat");
   	strcpy(t2_filename,"diff_times.dat");
   	strcpy(h0_filename,"heights0.dat");
   	strcpy(h1_filename,"heights1.dat");
   	strcpy(hits0_filename,"hits0.dat");
   	strcpy(hits1_filename,"hits1.dat");

   	FILE *ft0 = fopen(t0_filename, "w");
   	FILE *ft1 = fopen(t1_filename, "w");
   	FILE *ft2 = fopen(t2_filename, "w");
   	FILE *fh0 = fopen(h0_filename, "w");
   	FILE *fh1 = fopen(h1_filename, "w");
   	FILE *fhit0 = fopen(hits0_filename, "w");
   	FILE *fhit1 = fopen(hits1_filename, "w");
*/
   	// open the binary waveform file
   	FILE *f = fopen(filename, "r");
   	if (f == NULL) {
      		printf("Cannot find file \'%s\'\n", filename);
      		return 0;
   	}

   // read time header
   fread(&th, sizeof(th), 1, f);
//   printf("Found data for board #%d\n", th.board_serial_number);

   // read time bin widths
   memset(bin_width, sizeof(bin_width), 0);
   for (ch=0 ; ch<5 ; ch++) {
      fread(hdr, sizeof(hdr), 1, f);
      if (hdr[0] != 'C') {
         // event header found
         fseek(f, -4, SEEK_CUR);
         break;
      }
      i = hdr[3] - '0' - 1;
  //    printf("Found timing calibration for channel #%d\n", i+1);
      fread(&bin_width[i][0], sizeof(float), 1024, f);
   }
   
   // initialize statistics
   ndt = 0;
   sumdt = sumdt2 = 0;
   
   // Read to the first requested event
	int start_e=0,end_e=1000000000; 
	
	if (vm.count("start-event")) start_e = vm["start-event"].as<int>();
	if (vm.count("end-event")) end_e = vm["end-event"].as<int>();

	if (end_e <= start_e){
		std::cout << "Error : Requested start event index is after the end index" << std::endl;
		std::cout << "Error : Start event index is " << start_e << std::endl;
		std::cout << "Error : End event index is " << end_e << std::endl;
		return 0;		
	}




	// Main loop 
	drs4_run r(th);

	for (n= 0 ; n<start_e; n++) {
		i = fread(&eh, sizeof(eh), 1, f);
      		if (i < 1){
			std::cout << "Error : Requested start event index is beyond the end of the file" << std::endl;
			return 0;			
		}
		for (ch=0 ; ch<5 ; ch++) {
         		i = fread(hdr, sizeof(hdr), 1, f);
         		if (i < 1) break;
         		if (hdr[0] != 'C') {
            			// event header found
            			fseek(f, -4, SEEK_CUR);
            			break;
         		}
         		chn_index = hdr[3] - '0' - 1;
         		fread(voltage, sizeof(short), 1024, f);
		}   
	}




   // loop over all events in the data file
   for (n=start_e;n<end_e;n++) {

      // read event header
      i = fread(&eh, sizeof(eh), 1, f);
      if (i < 1) break;
      
	if((n%1000)==10){
		printf("# Event #%d\n", eh.event_serial_number);

		if(vm["run-number"].as<bool>())	{
			std::cout << "# Event Number : "; 
			r.print_event_number();
		}

		if(vm["run-rate"].as<bool>()){
			std::cout << "# Event Rate   : ";
			r.print_event_rate();
		}
	}

	// Save into event class structure
	drs4_event ev(eh,n,lower_threshold,upper_threshold);

      
      // reach channel data
      for (ch=0 ; ch<5 ; ch++) {
         i = fread(hdr, sizeof(hdr), 1, f);
         if (i < 1)
            break;
         if (hdr[0] != 'C') {
            // event header found
            fseek(f, -4, SEEK_CUR);
            break;
         }
         chn_index = hdr[3] - '0' - 1;
         fread(voltage, sizeof(short), 1024, f);
         
         for (i=0 ; i<1024 ; i++) {
            // convert data to volts
            waveform[chn_index][i] = (voltage[i] / 65536. - 0.5);
            
            // calculate time for this cell
            for (j=0,time[chn_index][i]=0 ; j<i ; j++)
               time[chn_index][i] += bin_width[chn_index][(j+eh.trigger_cell) % 1024];
         }
      }
      
      
      // align cell #0 of all channels
      t1 = time[0][(1024-eh.trigger_cell) % 1024];
      for (ch=1 ; ch<4 ; ch++) {
         t2 = time[ch][(1024-eh.trigger_cell) % 1024];
         dt = t1 - t2;
         for (i=0 ; i<1024 ; i++)
            time[ch][i] += dt;
      }
      
      t1 = t2 = 0;
//      threshold = 0.025;

	ev.fill(waveform,time);
	ev.fill_derived();

	r.add_event(ev);


	if( (vm["pulse-interevent"].as<bool>()) || (vm["pulse-hits"].as<bool>()) || (vm["pulse-times"].as<bool>()) || (vm["pulse-seq_times"].as<bool>()) || (vm["pulse-heights"].as<bool>())||(vm["pulse-delay"].as<bool>()) )
	{

	if(vm["pulse-interevent"].as<bool>()){
		if(channel>=0){
			ev.print_inter_times(channel);
		}else{
			ev.print_inter_times();	
		}
	}

	if(vm["pulse-delay"].as<bool>()){
		ev.print_delay_times();
	}

	if(vm["pulse-hits"].as<bool>()){
		if(channel>=0){
			ev.print_multiplicity(channel);
		}else{
			ev.print_multiplicity();	
		}
	}

	if(vm["pulse-times"].as<bool>()){
		if(channel>=0){
			ev.print_times(channel);
		}else{
			ev.print_times();	
		}
	}


	if(vm["pulse-seq_times"].as<bool>()){
			ev.print_seq_times();
	}

	if(vm["pulse-heights"].as<bool>()){
		if(channel>=0){
			ev.print_pulse_heights(channel);
		}else{
			ev.print_pulse_heights();	
		}
	}

	}else{

		//ev.print_event();
	}




      
      // find peak in channel 1 above threshold
/*      flag=0;
      nhits1=0;
      for (i=0 ; i<1022 ; i++){

         if (waveform[0][i] < threshold && waveform[0][i+1] >= threshold) {
            	t1 = (threshold-waveform[0][i])/(waveform[0][i+1]-waveform[0][i])*(time[0][i+1]-time[0][i])+time[0][i];
		nhits1++;
         }

         if (waveform[0][i] > pulse_height0) pulse_height0=waveform[0][i];

	 }

	 flag=0;
	 nhits2=0;

	// find peak in channel 2 above threshold
     for (i=0 ; i<1022 ; i++){

         if (waveform[1][i] < threshold && waveform[1][i+1] >= threshold) {
            	t2 = (threshold-waveform[1][i])/(waveform[1][i+1]-waveform[1][i])*(time[1][i+1]-time[1][i])+time[1][i];
		nhits2++;
         }

         if (waveform[1][i] > pulse_height1) pulse_height1=waveform[1][i];

	}
//	if(nhits1>0) printf("1 %d\n",nhits1);
//	if(nhits2>0) printf("2 %d\n",nhits2);

	if(nhits1>0){
	 	fprintf(ft0,"%1.3f\n",t1);
	 	fprintf(fh0,"%1.3f\n",pulse_height0);
	 	fprintf(fhit0,"%d\n",nhits1);
	}

	if(nhits2>0){
	 	fprintf(ft1,"%1.3f\n",t2);
	 	fprintf(fh1,"%1.3f\n",pulse_height1);
	 	fprintf(fhit1,"%d\n",nhits2);
	}

	if((nhits1>0) && (nhits2>0)){
	         ndt++;
	         dt = t2 - t1;
		 fprintf(ft2,"%1.3f\n",dt);
	}
*/


/*
    // calculate distance of peaks with statistics
      if (t1 > 0 || t2 > 0) {
         sumdt += dt;
         sumdt2 += dt*dt;

//	 printf("%1.3f %f %1.3f %f\n",t1,pulse_height0,t2,pulse_height1); 

	 fprintf(ft0,"%1.3f\n",t1);
	 fprintf(ft1,"%1.3f\n",t2);
	 fprintf(fh0,"%1.3f\n",pulse_height0);
	 fprintf(fh1,"%1.3f\n",pulse_height1);
	 fprintf(fhit0,"%d\n",nhits1);
	 fprintf(fhit1,"%d\n",nhits2);

      }
*/

//	pulse_height0=0;
//	pulse_height1=0;

   }
   
   // print statistics
   //printf("dT = %1.3lfns +- %1.1lfps\n", sumdt/ndt, 1000*sqrt(1.0/(ndt-1)*(sumdt2-1.0/ndt*sumdt*sumdt)));

	if(vm["run-times"].as<bool>())	r.print_event_times();
	if(vm["run-number"].as<bool>())	r.print_event_number();
	if(vm["run-rate"].as<bool>())	r.print_event_rate();
	r.print_event_timescale();
   
   return 1;
}

