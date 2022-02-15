#define SIZEOF_SINO 77*60
#define N_TX_BLK 120
#define N_TX_CRYS_PER_BLK 7
#define N_AX_CRYS_PER_BLK 12
#define N_AX_BLK 56 // 7 ax blocks per unit * 8 units
#define N_AX_Unit 84 //In Davis Listmode defination, we add virtual crystal ring for every 84 axial crystal, EX has 8 rings(672 crystal), we add 7 virtual crystal, so we have total 672+7=679 crystal
//#define N_AX_CRYS_PER_BLK 6
//#define N_AX_BLK 112 // 14 ax blocks per unit * 8 units
#define ADD_UNIT_GAP // comment out if no unit gap is desired
//#define Debug_mode // comment out if do not debug the code
#define N_TOF 27 //UIH defined [15,17,19,23,25,27], we use 27 to maintain compatibility with old format
#define N_TIME_BIN_PER_TOF 7
#include <vector>
#include <iostream>
//Fixed TOF binning BUG(-6 ~ +6 has been assigned into TOF bin 0 or middel bin) reported by Reimund
//Fixed TOF binning sign issue, zhxie@ucdavis.edu
//Fixed ADD_UNIT_GAP BUG(when add gap, N_AX_Unit should be 85 instead of 84) reported by Reimund

using namespace std;

struct Lut {
	short nv,nu;
};

struct Lm {
	short txIDA, axIDA, txIDB, axIDB, tof;
};

int main(int argc, char* argv[]) {

	// check arguments
	if (argc != 4) {
		cout << "Usage: " << argv[0] << " [fname_in] [fname_out] [fname_lut]" << endl;
		exit(1);
	}
        #ifdef ADD_UNIT_GAP
	cout << "ADD_UNIT_GAP is turned ON." << endl;
        #endif
        #ifdef Debug_mode
	cout << "Debug mode is turned ON." << endl;
        #endif

	cout << "Starting " << argv[0] << " ..." << endl;

	// read block sino lut
	string fname_lut = string(argv[3]);
	FILE* pfile_lut = fopen(fname_lut.c_str(), "rb");
	if (pfile_lut == NULL) {
		cout << fname_lut << " cannot be opened." << endl;
		exit(1);
	}
	Lut* plut = new Lut[SIZEOF_SINO]; // block sino lut
	cout << "Reading LUT "<< endl;
	fread(plut, sizeof(Lut), SIZEOF_SINO, pfile_lut);
	cout << "Done with LUT "<< endl;

	// set up block sino reverse lut
	vector< vector<int> > blk_idx(N_TX_BLK, vector<int>(N_TX_BLK, -1)); // initialize values to -1 to skip the Out-of-bound sinogram index
	for (int i = 0; i < SIZEOF_SINO; i++) {
	//cout << "Starting " <<plut[i].nv <<" ..."<<plut[i].nu<< " ..." <<i<< endl;  //Print transaxial ID (zhoujian's       'scanner.getDefaultSinogramCrystalPairs' format, remember to reverse TOF bin when swap LOR orientation)
		blk_idx[plut[i].nv-1][plut[i].nu-1] = i;
	//	blk_idx[plut[i].nu-1][plut[i].nv-1] = i;
	}
        cout << "BLK matrix size has  " << blk_idx.size() << "rows,and "<< blk_idx[0].size()<< "cols,\n each entry corresponds to blk# of LORA and blk#LORB in the transaxial dimension(range from 0~"<<N_TX_BLK-1<< ")" << endl;
	 
	cout << "Done with BLK index "<< endl;
	// declare input file
	string fname_in = string(argv[1]);
	FILE* pfile_in = fopen(fname_in.c_str(), "rb");
	if (pfile_in == NULL) {
		cout << fname_in << "cannot be opened." << endl;
		exit(1);
	}

	// declare output file
	string fname_out = string(argv[2]);
	FILE* pfile_out = fopen(fname_out.c_str(), "wb");
	if (pfile_out == NULL) {
		cout << fname_out << "cannot be opened." << endl;
		exit(1);
	}

	// main loop
        int read_count_outof_traxial_tof=1; //Because of zhoujian's  SinogramCrystalPairs has larger(22) accepted BLK Difference,which covers smaller FOV in transaxial plane, and TOF FOV is different
        int read_count_total=1; //
	Lm* plm = new Lm; // list-mode event
	unsigned long long* psino_blk = new unsigned long long[SIZEOF_SINO*N_AX_BLK*N_AX_BLK*N_TOF]; // block TOF sino
	while (!feof(pfile_in)) {

		int read_count = fread(plm, sizeof(Lm), 1, pfile_in); 

		for (int i = 0; i < read_count; i++) {
                        int TOF_AB;
			int txBiA = plm->txIDA / N_TX_CRYS_PER_BLK;
			int txBiB = plm->txIDB / N_TX_CRYS_PER_BLK;
			int TOF_AB_temp = plm->tof+3;
                        if (TOF_AB_temp>0||TOF_AB_temp==0){
			 TOF_AB = TOF_AB_temp/ N_TIME_BIN_PER_TOF;}
                        else{
                         TOF_AB = (TOF_AB_temp-N_TIME_BIN_PER_TOF)/ N_TIME_BIN_PER_TOF; }
                        // 
                       read_count_total++;
                       #ifdef ADD_UNIT_GAP
 			int axBiA_gap = plm->axIDA ;
			int axBiB_gap = plm->axIDB ;
			N_AX_Unit=85；
 			int axBiA = (axBiA_gap-axBiA_gap/N_AX_Unit)/ (N_AX_CRYS_PER_BLK);
			int axBiB = (axBiB_gap-axBiB_gap/N_AX_Unit)/ (N_AX_CRYS_PER_BLK);
                           #ifdef Debug_mode
                           if (i/1000000==1){
                           printf("Current LOR_A w ax_gap is: %d,LOR_B w ax_gap is %d,Tof bin is %d \n", axBiA_gap,axBiB_gap,plm->tof);
                           printf("Current LOR_A minus %d gap,LOR_B minus %d gap\n", axBiA_gap/N_AX_Unit,axBiB_gap/N_AX_Unit);
                           printf("Current LOR_A w/o ax_gap is: %d,LOR_B w/o ax_gap is %d\n", axBiA_gap-axBiA_gap/N_AX_Unit,axBiB_gap-axBiB_gap/N_AX_Unit);
                           printf("Current ABLK of LOR_A is: %d,ABLK of LOR_B  is %d, TOFBLK is %d \n \n", axBiA,axBiB,TOF_AB);}
                           #endif
                       #else //no gap insert in listmode
 			int axBiA = plm->axIDA / (N_AX_CRYS_PER_BLK);
			int axBiB = plm->axIDB / (N_AX_CRYS_PER_BLK);
                       #endif

                        

                        //cout << "Reading Listmode "<<txBiA<<"B "<<txBiB<< endl;
                        
			int idx_tx_blk = blk_idx[txBiA][txBiB]; // transaxial sinogram index ,120×120
			int idx_tx_blk_reverse = blk_idx[txBiB][txBiA]; // transaxial sinogram index ,120×120,reverse direction
                        int ind_blk_sino;
                        if (idx_tx_blk!=-1){
                         //The following part is to swap axial A and axial B for sinogram, this part is to match zhoujian's format(ZJ's forward projection). so axBiB is in front of axBiA
			ind_blk_sino = idx_tx_blk + SIZEOF_SINO * axBiB
				+ SIZEOF_SINO * N_AX_BLK * axBiA+ SIZEOF_SINO * N_AX_BLK * N_AX_BLK*(TOF_AB+13);} // 5-D sinogram 

                        else if(idx_tx_blk_reverse!=-1){   
                        // swap event A and event B's  position and reverse tof sign
                        TOF_AB=-TOF_AB;
                        idx_tx_blk=idx_tx_blk_reverse;
                        ind_blk_sino = idx_tx_blk + SIZEOF_SINO * axBiA
				+ SIZEOF_SINO * N_AX_BLK * axBiB+ SIZEOF_SINO * N_AX_BLK * N_AX_BLK*(TOF_AB+13);}

                        else  {                       
                        read_count_outof_traxial_tof++;

                         }
                       
                        //cout << "Reading ind_blk_sino "<<ind_blk_sino<< endl;
                        #ifdef Debug_mode
                        if (i/10==0){
                        if (ind_blk_sino<0||(TOF_AB+13)<0||(TOF_AB+13)>26||idx_tx_blk==(SIZEOF_SINO-1)||ind_blk_sino>(SIZEOF_SINO*N_AX_BLK*N_AX_BLK*N_TOF-1)){
	                cout << "Index Out-of-bound sinogram index detected. Skipping the following coincidence event: "<< endl;
                          printf("Current LOR crystal ID of tranxaxial plane is: %d,crystal ID of tranxaxial plane is %d,Tof bin is %d \n", plm->txIDA,plm->txIDB,plm->tof);
                          printf("Current BLK LORA tranxaxial ID is: %d,BLK LORB tranxaxial ID is %d(R:0~119),BLK Tof bin is %d  (R:0~26)\n", txBiA,txBiB,TOF_AB+13);
                          printf("Current BLK LORA Axial ID  is: %d, BLK LORB  Axial ID is %d \n", axBiA,axBiB);
                          printf("Current total index of LOR_A is: %d\n\n", ind_blk_sino);
                        }
                        }
                         
                        #endif
                        

			if (ind_blk_sino <0||ind_blk_sino>(SIZEOF_SINO*N_AX_BLK*N_AX_BLK*N_TOF-1)||(TOF_AB+13)<0||(TOF_AB+13)>26) { // provides fault tolerance
                       //ind_blk_sino <0 represents crystal-based UIH listmode out-of-bound zhoujian's sinogram,
                       //idx_tx_blk==SIZEOF_SINO*N_AX_BLK*N_AX_BLK*N_TOF-1 represents BLK-based listmodeID covered by zhoujian's sinogram, but some of crystal UIH listmode out-of-bound because of integer division in C is equivalent to matlab 'fix', which make the one last entry of sinogram has very high value .
                        read_count_outof_traxial_tof++;
			}
                        else {
			psino_blk[ind_blk_sino]++;
                        }
		}

	}

	printf("Totol event is %d, %d out of our transaxial or TOF FOV\n", read_count_total,read_count_outof_traxial_tof);
	// write to file
	cout << "Save sino... "<< endl;
	fwrite(psino_blk, sizeof(unsigned long long), SIZEOF_SINO*N_AX_BLK*N_AX_BLK*N_TOF, pfile_out);
	cout << "Done save sino." << endl;
	// cleanup

	delete[] plut;
	delete[] psino_blk;
	delete plm;
	fclose(pfile_lut);
	fclose(pfile_in);
  	fclose(pfile_out);
        // if segmentation fault, check each fclose and negative index

       
	return 0;
}
