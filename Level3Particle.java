import j4np.hipo5.data.Bank;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.lang.Math;

public class Level3Particle {

    float pid_resp=0;
    float track_resp=0;

    int PID=0;
    int MC_PID=0;

    int MC_PIndex=-1;
    int PIndex=-1;
    int Cal_Sector=0;
    int HTCC_Sector=0;
    int Sector=0;
    int Charge=0;

    Boolean show=false;

    int Track_index=-1;
    int[] track_clusterIDs = new int[6];
    float[] track_clusters = new float[6];
    double[] ecin_all_hits = new double[108];
    //double[] CF_out = new double[108];
    float[] CF_out = new float[]{0,0,0,0,0,0,0,0,0};

    //e-,pi+,pi-,e+,mu-,mu+
    int[] pid_label = new int[2]; //3

    double Px=0;
    double Py=0;
    double Pz=0;

    double P=0;
    double Theta=0;
    double Phi=0;

    double MC_Px=0;
    double MC_Py=0;
    double MC_Pz=0;

    double MC_P=0;
    double MC_Theta=0;
    double MC_Phi=0;

    int MC_Sector=0;

    double Nphe=0;
    float[] HTCC_adcs = new float[8];

    List<Integer> Cal_index=new ArrayList<Integer>();
    //double Cal_index=0;

    double SF=0;
    double ECAL_energy=0;

    double PCAL_energy=0;
    double ECIN_energy=0;
    double ECOUT_energy=0;

    double PCALLU=0;
    double PCALLV=0;
    double PCALLW=0;

    double ECINLU=0;
    double ECINLV=0;
    double ECINLW=0;

    double ECOUTLU=0;
    double ECOUTLV=0;
    double ECOUTLW=0;

    double ECOUTM2U=0;
    double ECOUTM2V=0;
    double ECOUTM2W=0;

    double PCALM2U=0;
    double PCALM2V=0;
    double PCALM2W=0;

    double ECINM2U=0;
    double ECINM2V=0;
    double ECINM2W=0;

    //below are values found in ECAL::adc bank
    //using cluster finder
    //energy is in ADC
    float PCAL_energy_fcf=0;
    float ECIN_energy_fcf=0;
    float ECOUT_energy_fcf=0;

    float PCALLU_fcf=-2;
    float PCALLV_fcf=-2;
    float PCALLW_fcf=-2;

    float ECINLU_fcf=-2;
    float ECINLV_fcf=-2;
    float ECINLW_fcf=-2;

    float ECOUTLU_fcf=-2;
    float ECOUTLV_fcf=-2;
    float ECOUTLW_fcf=-2;

    float PCALDU_fcf=0;
    float PCALDV_fcf=0;
    float PCALDW_fcf=0;

    float ECINDU_fcf=0;
    float ECINDV_fcf=0;
    float ECINDW_fcf=0;

    float ECOUTDU_fcf=0;
    float ECOUTDV_fcf=0;
    float ECOUTDW_fcf=0;

    double PCAL_scale=6.5;
    double ECAL_scale=8.0;

    double PCal_UMax_cut = 80.;
    double PCal_UMin_cut = 3.; // For outbendings this doesn't matter much;
    double PCal_VMax_cut = 73.;
    double PCal_VMin_cut = 3.;
    double PCal_WMax_cut = 73.;
    double PCal_WMin_cut = 3.;

    double ECin_UMax_cut = 35.;
    double ECin_UMin_cut = 2.5; // For outbendings this doesn't matter much;
    double ECin_VMax_cut = 34.;
    double ECin_VMin_cut = 3.;
    double ECin_WMax_cut = 34.;
    double ECin_WMin_cut = 3.;

    double ECout_UMax_cut = 34.;
    double ECout_UMin_cut = 2; // For outbendings this doesn't matter much;
    double ECout_VMax_cut = 34;
    double ECout_VMin_cut = 3.;
    double ECout_WMax_cut = 34.;
    double ECout_WMin_cut = 3.;

    public Level3Particle(){

    }

    public void setTrackResp(float r){
        track_resp=r;
    }

    public void setPidResp(float r){
        pid_resp=r;
    }

    public void setShow(Boolean sh){
        show=sh;
    }

    //using ECAL::adc
    public void read_cal_bank_from_cf_pred(float[] cf_pred, int[] cf_strips, Bank ECAL_Bank){
        PCALLU_fcf=Math.abs(cf_pred[0]);
        PCALLV_fcf=Math.abs(cf_pred[1]);
        PCALLW_fcf=Math.abs(cf_pred[2]);
        ECINLU_fcf=Math.abs(cf_pred[3]);
        ECINLV_fcf=Math.abs(cf_pred[4]);
        ECINLW_fcf=Math.abs(cf_pred[5]);
        ECOUTLU_fcf=Math.abs(cf_pred[6]);
        ECOUTLV_fcf=Math.abs(cf_pred[7]);
        ECOUTLW_fcf=Math.abs(cf_pred[8]);

        for(int k = 0; k < ECAL_Bank.getRows(); k++){
            
            int   sect = ECAL_Bank.getInt("sector", k);
            int  layer = ECAL_Bank.getInt("layer", k);
            int  strip = ECAL_Bank.getInt("component", k);
            int    ADC = ECAL_Bank.getInt("ADC", k);

            if(ADC>0.0){
                
                //----------------
                if(sect==Sector){
                    if(layer==1){
                        if(strip>(cf_strips[0]-6) && strip<(cf_strips[0]+6)){
                            PCALDU_fcf++;
                            PCAL_energy_fcf+=ADC;
                        }
                    } else if(layer==2){
                        if(strip>(cf_strips[1]-6) && strip<(cf_strips[1]+6)){
                            PCALDV_fcf++;
                            PCAL_energy_fcf+=ADC;
                        }
                    } else if(layer==3){
                        if(strip>(cf_strips[2]-6) && strip<(cf_strips[2]+6)){
                            PCALDW_fcf++;
                            PCAL_energy_fcf+=ADC;
                        }
                    } else if(layer==4){
                        if(strip>(cf_strips[3]-6) && strip<(cf_strips[3]+6)){
                            ECINDU_fcf++;
                            ECIN_energy_fcf+=ADC;
                        }
                    } else if(layer==5){
                        if(strip>(cf_strips[4]-6) && strip<(cf_strips[4]+6)){
                            ECINDV_fcf++;
                            ECIN_energy_fcf+=ADC;
                        }
                    } else if(layer==6){
                        if(strip>(cf_strips[5]-6) && strip<(cf_strips[5]+6)){
                            ECINDW_fcf++;
                            ECIN_energy_fcf+=ADC;
                        }
                    } else if(layer==7){
                        if(strip>(cf_strips[6]-6) && strip<(cf_strips[6]+6)){
                            ECOUTDU_fcf++;
                            ECOUT_energy_fcf+=ADC;
                        }
                    } else if(layer==8){
                        if(strip>(cf_strips[7]-6) && strip<(cf_strips[7]+6)){
                            ECOUTDV_fcf++;
                            ECOUT_energy_fcf+=ADC;
                        }
                    } else if(layer==9){
                        if(strip>(cf_strips[8]-6) && strip<(cf_strips[8]+6)){
                            ECOUTDW_fcf++;
                            ECOUT_energy_fcf+=ADC;
                        }
                    }
                }
            }
        }
        

    }

    //using REC::Track
    public void find_sector_track(Bank TrackBank){
        if(show){
            System.out.println("REC::Track");
            TrackBank.show();
        }
        
        for (int k = 0; k < TrackBank.getRows(); k++) {
            int pindex = TrackBank.getInt("pindex", k);
            int ind = TrackBank.getInt("index", k);
            int sectorTrk = TrackBank.getInt("sector", k);
            if(pindex==PIndex){
                Sector=sectorTrk;
                Track_index=ind;
            }
        }
    }

    //using HTCC:adc
    public void find_HTCC_ADCs(Bank HTCCBank){
        for(int k = 0; k < HTCCBank.getRows(); k++){
            int   sect = HTCCBank.getInt("sector", k);
            int  layer = HTCCBank.getInt("layer", k); //1 or 2
            int  component = HTCCBank.getInt("component", k); //1-4
            int    ADC = HTCCBank.getInt("ADC", k);

            //double energy = (ADC / 5000.0); // &&energy<1.0

            if (ADC >= 0.0 && sect == Sector) {
                int index = ((layer - 1) * 4 + component) - 1;
                if(index>=0 && index<8){
                    HTCC_adcs[index]=ADC;
                }   
            }
        }
    }

    //using TimeBasedTrkg::TBTracks or HitBasedTrkg::HBTracks
    public void find_track_clusterIDs(Bank TrackBank){
        if(show){
            System.out.println("HB/TB::HB/TBTrack");
            System.out.printf("Track index %d \n",Track_index);
            TrackBank.show();
        }
        //need to already find track sector from REC::Track
        int sect = TrackBank.getInt("sector", Track_index);
        if(sect==Sector){
            for(int i=0;i<6;i++){
                track_clusterIDs[i]=TrackBank.getInt("Cluster"+String.valueOf(i+1)+"_ID", Track_index);
            }
        }
    }

    //using TimeBasedTrkg::TBClusters or HitBasedTrkg::HBClusters
    public void find_track_clusters(Bank TrackBank){
        if(show){
            System.out.println("HB/TB::HB/TBTrack");
            System.out.println("Recorded Cluster IDs"+Arrays.toString(track_clusterIDs));
            TrackBank.show();
        }

        //need to already find track sector from REC::Track
        for (int k = 0; k < TrackBank.getRows(); k++) {
            int id = TrackBank.getInt("id", k);
            for(int i=0;i<6;i++){
                if(track_clusterIDs[i]==id){
                    track_clusters[i]=TrackBank.getFloat("avgWire", k)/112;
                }
            }
        }
    }

    public void find_sector_cal(Bank CalBank){
        for (int k = 0; k < CalBank.getRows(); k++) {
            int pindex = CalBank.getInt("pindex", k);
            int sectorTrk = CalBank.getInt("sector", k);
            if(pindex==PIndex){Cal_Sector=sectorTrk;}
        }
    }

    public void read_HTCC_bank(Bank HTCCBank){
        for (int k = 0; k < HTCCBank.getRows(); k++) {
            int pindex = HTCCBank.getInt("pindex", k);
            double nphe = HTCCBank.getFloat("nphe", k);
            int sector = HTCCBank.getInt("sector", k);
            if(pindex==PIndex){
                Nphe=nphe;
                HTCC_Sector=sector;
            }
        }
    }

    public Boolean isMip(Bank ECAL_Bank){
        Boolean ismip=true;
        for (int k = 0; k < ECAL_Bank.getRows(); k++) {
            float widthU = ECAL_Bank.getFloat("widthU",k);
            float widthV = ECAL_Bank.getFloat("widthV",k);
            float widthW = ECAL_Bank.getFloat("widthW",k);
            if(Cal_index.contains(k)){ 
            //if(Cal_index==k){
                if(widthU>1 && widthV>1 && widthW>1){
                    ismip=false;
                }
            }

        }
        return ismip;
    }

    //using ECAL clusters
    public Boolean check_FID_Cal_Clusters(Bank ECAL_Bank){
        Boolean Fid=true;
        //ECAL_Bank.show();
        //System.out.printf("index %d sector %d\n",PIndex,Sector);
        for (int k = 0; k < ECAL_Bank.getRows(); k++) {
            int sector=ECAL_Bank.getInt("sector", k);
            int layer=ECAL_Bank.getInt("layer",k);
            double u=(double) ECAL_Bank.getInt("coordU",k);
            double v=(double) ECAL_Bank.getInt("coordV",k);
            double w=(double) ECAL_Bank.getInt("coordW",k);
            double energy=ECAL_Bank.getFloat("energy",k);
            
            if(layer==1 && Cal_index.contains(k)){// Cal_index==k){
                u=u/PCAL_scale;
                v=v/PCAL_scale;
                w=w/PCAL_scale;
                /*System.out.printf("layer %d, sector %d, u %f,v %f,w %f,energy %f\n",layer,sector,u,v,w,energy);
                System.out.printf("Sector %d, energy PCAL %f, ECIN %f , ECOUT %f\n",Sector,PCAL_energy,ECIN_energy,ECOUT_energy);*/
                if(u>PCal_UMax_cut||u<PCal_UMin_cut){Fid=false;}
                if(v>PCal_VMax_cut||v<PCal_VMin_cut){Fid=false;}
                if(w>PCal_WMax_cut||w<PCal_WMin_cut){Fid=false;}
            } else if(layer==4 && Cal_index.contains(k)){// Cal_index==k){
                u=u/ECAL_scale;
                v=v/ECAL_scale;
                w=w/ECAL_scale;
                /*System.out.printf("layer %d, sector %d, u %f,v %f,w %f,energy %f\n",layer,sector,u,v,w,energy);
                System.out.printf("Sector %d, energy PCAL %f, ECIN %f , ECOUT %f\n",Sector,PCAL_energy,ECIN_energy,ECOUT_energy);*/
                if(u>ECin_UMax_cut||u<ECin_UMin_cut){Fid=false;}
                if(v>ECin_VMax_cut||v<ECin_VMin_cut){Fid=false;}
                if(w>ECin_WMax_cut||w<ECin_WMin_cut){Fid=false;}
            } else if(layer==7 &&Cal_index.contains(k)){// Cal_index==k){
                u=u/ECAL_scale;
                v=v/ECAL_scale;
                w=w/ECAL_scale;
                /*System.out.printf("layer %d, sector %d, u %f,v %f,w %f,energy %f\n",layer,sector,u,v,w,energy);
                System.out.printf("Sector %d, energy PCAL %f, ECIN %f , ECOUT %f\n",Sector,PCAL_energy,ECIN_energy,ECOUT_energy);*/
                if(u>ECout_UMax_cut||u<ECout_UMin_cut){Fid=false;}
                if(v>ECout_VMax_cut||v<ECout_VMin_cut){Fid=false;}
                if(w>ECout_WMax_cut||w<ECout_WMin_cut){Fid=false;}
            }
            
        }
        return Fid;
    }

    

    public Boolean check_Energy_Dep_Cut(){
        Boolean pass=false;
        //System.out.printf("energy P %f i %f o %f t %f\n",PCAL_energy,ECIN_energy,ECOUT_energy,tot_e);
        if(SF>0.2&&PCAL_energy>0.06){
            pass=true;
        }
        return pass;
    }

    public Boolean check_SF_cut(){
        Boolean pass=false;
        if(SF>0.2){pass=true;}
        return pass;
    }

    //using ECAL::adc
    public void fillECin(Bank ECAL_Bank){
        
        for(int k = 0; k < ECAL_Bank.getRows(); k++){
            
            int   sect = ECAL_Bank.getInt("sector", k);
            int  layer = ECAL_Bank.getInt("layer", k);
            int  strip = ECAL_Bank.getInt("component", k);
            int    ADC = ECAL_Bank.getInt("ADC", k);

            if(ADC>0.0){
                double energy = (ADC/10000.0)/1.5/3.0;
                //----------------
                if(sect==Sector){
                    if (layer > 3 && layer < 7) {
                        int index = ((layer-4)*36+strip)-1;

                        if (index < 108) {
                            ecin_all_hits[index]=energy;
                        }
                    }
                }
            }
        }
    }

    //REC::Calorimeter
    public void read_Cal_Bank(Bank ECAL_Bank){
        if(show){
            System.out.println("REC::Calorimeter");
            ECAL_Bank.show();
        }
        //ECAL_Bank.show();
        for (int k = 0; k < ECAL_Bank.getRows(); k++) {
            int pindex = ECAL_Bank.getInt("pindex", k);
            int index = ECAL_Bank.getInt("index", k);
            double energy = ECAL_Bank.getFloat("energy", k);
            int sector = ECAL_Bank.getInt("sector", k);
            int layer=ECAL_Bank.getInt("layer",k);
            float lu=ECAL_Bank.getFloat("lu",k);
            float lv=ECAL_Bank.getFloat("lv",k);
            float lw=ECAL_Bank.getFloat("lw",k);
            double m2u=ECAL_Bank.getFloat("m2u",k);
            double m2v=ECAL_Bank.getFloat("m2v",k);
            double m2w=ECAL_Bank.getFloat("m2w",k);
            if (pindex == PIndex) {
                Cal_index.add(index);
                //Cal_index=index;
                if(layer==1){
                    PCALLU=lu;
                    PCALLV=lv;
                    PCALLW=lw;
                    PCALM2U=m2u;
                    PCALM2V=m2v;
                    PCALM2W=m2w;
                    PCAL_energy=energy;
                    CF_out[0]=lu;
                    CF_out[1]=lv;
                    CF_out[2]=lw;
                } else if(layer==4){
                    ECINLU=lu;
                    ECINLV=lv;
                    ECINLW=lw;
                    ECINM2U=m2u;
                    ECINM2V=m2v;
                    ECINM2W=m2w;
                    ECIN_energy=energy;
                    CF_out[3]=lu;
                    CF_out[4]=lv;
                    CF_out[5]=lw;
                } else if(layer==7){
                    ECOUTLU=lu;
                    ECOUTLV=lv;
                    ECOUTLW=lw;
                    ECOUTM2U=m2u;
                    ECOUTM2V=m2v;
                    ECOUTM2W=m2w;
                    ECOUT_energy=energy;
                    CF_out[6]=lu;
                    CF_out[7]=lv;
                    CF_out[8]=lw;
                }
            }
        }
        ECAL_energy=PCAL_energy+ECIN_energy+ECOUT_energy;
        SF=ECAL_energy/P;


    }

    public Boolean TruthMatch(double Plim,double Thetalim, double Philim){
        Boolean truthmatched=false;
        double resP=Math.abs(P-MC_P);
        double resTheta=Math.abs(Theta-MC_Theta);
        double resPhi=Math.abs(Phi-MC_Phi);
        if(resP<Plim && resTheta<Thetalim && resPhi<Philim){
            truthmatched=true;
        }
        return truthmatched;

    }

    public void find_ClosestMCParticle(Bank MCPartBank){
        double min_resPx=9999;
        double min_resPy=9999;
        double min_resPz=9999;
        int best_ind=-1;
        for (int i = 0; i < MCPartBank.getRows(); i++) {
            read_MCParticle_Bank(i, MCPartBank);
            double resPx=Math.abs(Px-MC_Px);
            double resPy=Math.abs(Py-MC_Py);
            double resPz=Math.abs(Pz-MC_Pz);
            if(resPx<min_resPx && resPy<min_resPy && resPz<min_resPz){
                min_resPx=resPx;
                min_resPy=resPy;
                min_resPz=resPz;
                best_ind=i;
            }
        }
        read_MCParticle_Bank(best_ind, MCPartBank);
        MC_PIndex=best_ind;
    }

    public void read_MCParticle_Bank(int pindex, Bank PartBank) {
        MC_PID = PartBank.getInt("pid", pindex);
        MC_Px = PartBank.getFloat("px", pindex);
        MC_Py = PartBank.getFloat("py", pindex);
        MC_Pz = PartBank.getFloat("pz", pindex);
        MC_P = Math.sqrt(MC_Px * MC_Px + MC_Py * MC_Py + MC_Pz * MC_Pz);
        MC_Theta = Math.acos(MC_Pz / MC_P);// Math.atan2(Math.sqrt(px*px+py*py),pz);
        MC_Phi = Math.atan2(MC_Py, MC_Px);
        double theta_deg=MC_Theta* (180.0 / Math.PI);
        double phi_deg=MC_Phi* (180.0 / Math.PI);
        if(theta_deg>5 && theta_deg<35){
            if(MC_PID==22 || MC_PID==2112){
                if (phi_deg > -30 && phi_deg < 30) {
                    MC_Sector = 1;
                } else if (phi_deg > 30 && phi_deg < 90) {
                    MC_Sector = 2;
                } else if (phi_deg > 90 && phi_deg < 150) {
                    MC_Sector = 3;
                } else if (phi_deg > 150 && phi_deg < -150) {
                    MC_Sector = 4;
                } else if (phi_deg > -150 && phi_deg < -90) {
                    MC_Sector = 5;
                } else if (phi_deg > -90 && phi_deg < -30) {
                    MC_Sector = 6;
                }
            } else{
                //shoudl this be different for charged particles?
                // they'll drift in phi
                //this also depends if negative or positive charge
                //and depends on field...
                //done below for e- in inbending
                if (phi_deg > -15 && phi_deg < 40) {
                    MC_Sector = 1;
                } else if (phi_deg > 45 && phi_deg < 90) {
                    MC_Sector = 2;
                } else if (phi_deg > 105 && phi_deg < 145) {
                    MC_Sector = 3;
                } else if (phi_deg > 165 && phi_deg < -150) {
                    MC_Sector = 4;
                } else if (phi_deg > -135 && phi_deg < -90) {
                    MC_Sector = 5;
                } else if (phi_deg > -75 && phi_deg < -30) {
                    MC_Sector = 6;
                }

                //trying to account for regions of no acceptance...
                /*if (phi_deg > -20 && phi_deg < 40) {
                    MC_Sector = 1;
                } else if (phi_deg > 40 && phi_deg < 100) {
                    MC_Sector = 2;
                } else if (phi_deg > 100 && phi_deg < 160) {
                    MC_Sector = 3;
                } else if (phi_deg > 160 && phi_deg < -140) {
                    MC_Sector = 4;
                } else if (phi_deg > -140 && phi_deg < -80) {
                    MC_Sector = 5;
                } else if (phi_deg > -80 && phi_deg < -20) {
                    MC_Sector = 6;
                } */
                
            }
        }
        MC_PIndex=pindex;
        
    }

    public void find_ClosestRECParticle(Bank PartBank){
        double min_resPx=9999;
        double min_resPy=9999;
        double min_resPz=9999;
        int best_ind=-1;
        for (int i = 0; i < PartBank.getRows(); i++) {
            read_Particle_Bank(i, PartBank);
            double resPx=Math.abs(Px-MC_Px);
            double resPy=Math.abs(Py-MC_Py);
            double resPz=Math.abs(Pz-MC_Pz);
            if(resPx<min_resPx && resPy<min_resPy && resPz<min_resPz){
                min_resPx=resPx;
                min_resPy=resPy;
                min_resPz=resPz;
                best_ind=i;
            }
        }
        read_Particle_Bank(best_ind, PartBank);
        PIndex=best_ind;
    }

    public void read_Particle_Bank(int pindex, Bank PartBank) {
        if(show){
            System.out.printf("\n New Particle pindex %d \n",pindex);
            System.out.println("REC::Particle");
            PartBank.show();
        }
        int pid = PartBank.getInt("pid", pindex);
        int status = PartBank.getInt("status", pindex);
        int charge = PartBank.getInt("charge", pindex);
        double px = PartBank.getFloat("px", pindex);
        double py = PartBank.getFloat("py", pindex);
        double pz = PartBank.getFloat("pz", pindex);
        if (Math.abs(status) >= 2000 && Math.abs(status) < 4000) {
            PID = pid;
            Charge = charge;
            Px = px;
            Py = py;
            Pz = pz;
            P = Math.sqrt(px * px + py * py + pz * pz);
            Theta = Math.acos(pz / P);// Math.atan2(Math.sqrt(px*px+py*py),pz);
            Phi = Math.atan2(py, px);
            PIndex = pindex;
        }
    }

    //e-,e+,pi-,pi+,mu-,mu+
    public int set_pid_label(){

        //System.out.printf("MC PID %d ",MC_PID);

        int out_ind=-1;
        if(MC_PID==11){
            out_ind=0;
        } else if(MC_PID==-211){
            out_ind=1;
        }/*else if(MC_PID==211){
            out_ind=2;
        }else if(MC_PID==-11){
            out_ind=3;
        }else if(MC_PID==-13){
            out_ind=4;
        }else if(MC_PID==13){
            out_ind=5;
        }*/

        //System.out.printf("out_ind %d ",out_ind);

        for(int i=0;i<2;i++){
            if(i==out_ind){pid_label[i]=1;}
            else{pid_label[i]=0;}
        }

        return out_ind;
    }

    //35 floating-point values and 6 integer values...
    public String get_csv_out(){
        StringBuilder csvLineBuilder = new StringBuilder();

        csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,", PCAL_energy/2.0, ECIN_energy/2.0, ECOUT_energy/2.0));// /3
        csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,", PCALLU/500.0, PCALLV/500.0, PCALLW/500.0)); // /2000
        csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,", ECINLU/500.0, ECINLV/500.0, ECINLW/500.0));
        csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,", ECOUTLU/500.0, ECOUTLV/500.0, ECOUTLW/500.0));
        csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,", PCALM2U/1000.0, PCALM2V/1000.0, PCALM2W/1000.0)); // 
        csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,", ECINM2U/1000.0, ECINM2V/1000.0, ECINM2W/1000.0));
        csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,", ECOUTM2U/1000.0, ECOUTM2V/1000.0, ECOUTM2W/1000.0));
        
        for (float value : track_clusters) {
            csvLineBuilder.append(String.format("%.6f,", value));
        }
        
        for (double value : HTCC_adcs) {
            csvLineBuilder.append(String.format("%.6f,",value/30000));
        }
        
        for (int value : pid_label) {
            csvLineBuilder.append(String.format("%d,", value));
        }

        csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,", P, Theta*(180.0 / Math.PI), Phi*(180.0 / Math.PI)));

        // Remove the trailing comma
        csvLineBuilder.deleteCharAt(csvLineBuilder.length() - 1);
        // Convert StringBuilder to String
        return csvLineBuilder.toString();
    }

    //35 floating-point values and 6 integer values...
    public String get_csv_out_fromCFPred(){
        StringBuilder csvLineBuilder = new StringBuilder();

        csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,", PCAL_energy_fcf/150000.0, ECIN_energy_fcf/150000.0, ECOUT_energy_fcf/150000.0));// /3
        csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,", PCALLU_fcf/500.0, PCALLV_fcf/500.0, PCALLW_fcf/500.0)); // /2000
        csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,", ECINLU_fcf/500.0, ECINLV_fcf/500.0, ECINLW_fcf/500.0));
        csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,", ECOUTLU_fcf/500.0, ECOUTLV_fcf/500.0, ECOUTLW_fcf/500.0));
        csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,", PCALDU_fcf/16.0, PCALDV_fcf/16.0, PCALDW_fcf/16.0)); //loop over at most 11 strips
        csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,", ECINDU_fcf/16.0, ECINDV_fcf/16.0, ECINDW_fcf/16.0)); //5 each side of pred
        csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,", ECOUTDU_fcf/16.0, ECOUTDV_fcf/16.0, ECOUTDW_fcf/16.0));
        
        for (float value : track_clusters) {
            csvLineBuilder.append(String.format("%.6f,", value));
        }
        
        for (double value : HTCC_adcs) {
            csvLineBuilder.append(String.format("%.6f,",value/35000));
        }
        
        for (int value : pid_label) {
            csvLineBuilder.append(String.format("%d,", value));
        }

        csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,", P, Theta*(180.0 / Math.PI), Phi*(180.0 / Math.PI)));

        // Remove the trailing comma
        csvLineBuilder.deleteCharAt(csvLineBuilder.length() - 1);
        // Convert StringBuilder to String
        return csvLineBuilder.toString();
    }

    //35 floating-point values and 6 integer values...
    public float[] get_vars_forpid(){
        float[] vars_for_pid = new float[35];
        vars_for_pid[0]=PCAL_energy_fcf/150000;
        vars_for_pid[1]=ECIN_energy_fcf/150000;
        vars_for_pid[2]=ECOUT_energy_fcf/150000;
        vars_for_pid[3]=PCALLU_fcf/500;
        vars_for_pid[4]=PCALLV_fcf/500;
        vars_for_pid[5]=PCALLW_fcf/500;
        vars_for_pid[6]=ECINLU_fcf/500;
        vars_for_pid[7]=ECINLV_fcf/500;
        vars_for_pid[8]=ECINLW_fcf/500;
        vars_for_pid[9]=ECOUTLU_fcf/500;
        vars_for_pid[10]=ECOUTLV_fcf/500;
        vars_for_pid[11]=ECOUTLW_fcf/500;
        vars_for_pid[12]=PCALDU_fcf/16;
        vars_for_pid[13]=PCALDV_fcf/16;
        vars_for_pid[14]=PCALDW_fcf/16;
        vars_for_pid[15]=ECINDU_fcf/16;
        vars_for_pid[16]=ECINDV_fcf/16;
        vars_for_pid[17]=ECINDW_fcf/16;
        vars_for_pid[18]=ECOUTDU_fcf/16;
        vars_for_pid[19]=ECOUTDV_fcf/16;
        vars_for_pid[20]=ECOUTDW_fcf/16;
        
        int n=21;
        for (float value : track_clusters) {
            vars_for_pid[n]=value;
            n++;
        }
        
        for (float value : HTCC_adcs) {
            vars_for_pid[n]=value/35000;
            n++;
        }
        
        return vars_for_pid;
    }

    public String get_csv_cf_out(){
        StringBuilder csvLineBuilder = new StringBuilder();
        // Append values from track_clusters array to the StringBuilder
        for (float value : track_clusters) {
            csvLineBuilder.append(String.format("%.6f,", value));
        }

        // Append values from CF_out array to the StringBuilder
        for (double value : CF_out) {
            csvLineBuilder.append(String.format("%.6f,", value/500)); //500 using LU/LV/LW
        }
        csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,", P, Theta*(180.0 / Math.PI), Phi*(180.0 / Math.PI)));
        // Remove the trailing comma
        csvLineBuilder.deleteCharAt(csvLineBuilder.length() - 1);
        // Convert StringBuilder to String
        return csvLineBuilder.toString();
    }

    public double getM(int pid){
        switch(pid){
            case 22: return 0;
            case 11: return 0.000511;
            case -11: return 0.000511;
            case 211: return 0.13957;
            case -211: return 0.13957;
            case 13: return 0.10566;
            case -13: return 0.10566;
            case 321: return 0.49368;
            case -321: return 0.49368;
            case 2212: return 0.938272;
            case 2112: return 0.939565;
            default: return -1;
        }
    }

}
