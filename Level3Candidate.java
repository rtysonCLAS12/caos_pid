import j4np.hipo5.data.Bank;
import j4np.hipo5.data.CompositeNode;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.lang.Math;

import j4np.neural.classifier.NeuralClassifierModel;

public class Level3Candidate {

    float pid_resp=0;
    float track_resp=0;

    int PID=-1;

    int PIndex=-1;
    int Cal_Sector=0;
    int Sector=0;
    int Charge=0;
    int Pred_Charge=0;

    Boolean show=false;
    Boolean unmatched=true;

    float[] track_clusters = new float[6];
    double[] ecin_all_hits = new double[108];
    //double[] CF_out = new double[108];
    float[] CF_out = new float[]{0,0,0,0,0,0,0,0,0,0,0};

    //e-,pi+,pi-,e+,mu-,mu+
    int[] pid_label = new int[2]; //3

    double Vx=0;
    double Vy=0;
    double Vz=0;
    double Px=0;
    double Py=0;
    double Pz=0;

    double chi2pid=0;

    double P=0;
    double Theta=0;
    double Phi=0;

    double Pred_Px=0;
    double Pred_Py=0;
    double Pred_Pz=0;
    float Pred_P=0;
    float Pred_Theta=0;
    float Pred_Phi=0;

    double Nphe=0;
    float[] HTCC_adcs = new float[]{0,0,0,0,0,0,0,0};

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

    double FTOF_comp=0;
    double FTOF_path=0;

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

    public Level3Candidate(){

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

    public void print(){
        System.out.println("");
        System.out.printf("PID %d pindex %d sector %d cal_sector %d charge %d\n",PID,PIndex,Sector,Cal_Sector,Charge);
        System.out.printf("pred Px %f Py %f Pz %f \n",Pred_Px,Pred_Py,Pred_Pz);
        System.out.println("Clusters ");
        for(int k = 0; k < 6; k++) System.out.printf(" %f ", track_clusters[k]);
        System.out.println("");
        System.out.printf("PCAL E %f ECIN # %f ECOUT E %f\n",PCAL_energy,ECIN_energy,ECOUT_energy);
        System.out.println("ECAL positions ");
        for(int k = 0; k < 9; k++) System.out.printf(" %f ", CF_out[k]);
        System.out.println("");
        System.out.println("HTCC ADCs ");
        for(int k = 0; k < 8; k++) System.out.printf(" %f ", HTCC_adcs[k]);
        System.out.println("");


    }
    

    public void printbanks(CompositeNode tr,CompositeNode pt, Bank RECb, Bank CALb, Bank HTCCb){
        for(int i=0;i<tr.getRows();i++){
            System.out.printf("\n New Track index %d \n",i);
            System.out.printf("sector = %2d, charge = %3d segments [ ", tr.getInt(1,i),tr.getInt(2,i));
            for(int k = 0; k < 6; k++) System.out.printf(" %9.5f ",tr.getDouble(10+k,i));
            System.out.printf("] p %9.5f  %9.5f %9.5f\n",pt.getDouble(4,i),pt.getDouble(5,i),pt.getDouble(6,i)); 
        }

        System.out.println("\nREC::Particle");
        RECb.show();
        System.out.println("\nREC::Calorimeter");
        CALb.show();
        System.out.println("\nHTCC ADC");
        HTCCb.show();
        

    }

    //using composite nodes 3210,2 and 3210,3
    public void find_AI_tracks(CompositeNode tr,CompositeNode pt, int i){
        if(show){
            System.out.printf("\n New Track index %d \n",i);
            System.out.printf("sector = %2d, charge = %3d segments [ ", tr.getInt(1,i),tr.getInt(2,i));
            for(int k = 0; k < 6; k++) System.out.printf(" %9.5f ",tr.getDouble(10+k,i));
            System.out.printf("] p %9.5f  %9.5f %9.5f\n",pt.getDouble(4,i),pt.getDouble(5,i),pt.getDouble(6,i));
        }

        Sector=tr.getInt(1,i);
        Pred_Charge=tr.getInt(2,i);
        for(int k = 0; k < 6; k++){
            track_clusters[k]=(float) tr.getDouble(10+k,i)/112;
        }
        Pred_Px=pt.getDouble(4,i);
        Pred_Py=pt.getDouble(5,i);
        Pred_Pz=pt.getDouble(6,i);
        Pred_P = (float) Math.sqrt(Pred_Px * Pred_Px + Pred_Py * Pred_Py + Pred_Pz * Pred_Pz);
        Pred_Theta = (float) Math.acos((float)Pred_Pz / Pred_P);// Math.atan2(Math.sqrt(Pred_Px*Pred_Px+Pred_Py*Pred_Py),Pred_Pz);
        Pred_Phi = (float) Math.atan2(Pred_Py, Pred_Px);

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
                        if(strip>(cf_strips[0]-4) && strip<(cf_strips[0]+4)){
                            PCALDU_fcf++;
                            PCAL_energy_fcf+=ADC;
                        }
                    } else if(layer==2){
                        if(strip>(cf_strips[1]-4) && strip<(cf_strips[1]+4)){
                            PCALDV_fcf++;
                            PCAL_energy_fcf+=ADC;
                        }
                    } else if(layer==3){
                        if(strip>(cf_strips[2]-4) && strip<(cf_strips[2]+4)){
                            PCALDW_fcf++;
                            PCAL_energy_fcf+=ADC;
                        }
                    } else if(layer==4){
                        if(strip>(cf_strips[3]-4) && strip<(cf_strips[3]+4)){
                            ECINDU_fcf++;
                            ECIN_energy_fcf+=ADC;
                        }
                    } else if(layer==5){
                        if(strip>(cf_strips[4]-4) && strip<(cf_strips[4]+4)){
                            ECINDV_fcf++;
                            ECIN_energy_fcf+=ADC;
                        }
                    } else if(layer==6){
                        if(strip>(cf_strips[5]-4) && strip<(cf_strips[5]+4)){
                            ECINDW_fcf++;
                            ECIN_energy_fcf+=ADC;
                        }
                    } else if(layer==7){
                        if(strip>(cf_strips[6]-4) && strip<(cf_strips[6]+4)){
                            ECOUTDU_fcf++;
                            ECOUT_energy_fcf+=ADC;
                        }
                    } else if(layer==8){
                        if(strip>(cf_strips[7]-4) && strip<(cf_strips[7]+4)){
                            ECOUTDV_fcf++;
                            ECOUT_energy_fcf+=ADC;
                        }
                    } else if(layer==9){
                        if(strip>(cf_strips[8]-4) && strip<(cf_strips[8]+4)){
                            ECOUTDW_fcf++;
                            ECOUT_energy_fcf+=ADC;
                        }
                    }
                }
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

    public void applyTriangCut(){
        double SFPCAL=PCAL_energy/P;
        double SFECIN=ECIN_energy/P;
        double SFECOUT=ECOUT_energy/P;
        
        if(PID==11 && P>4.5){
            if(SFECIN > (0.2 - SFPCAL)){
                //all good
            } else{
                PID=-211;
            }
        }
    }

    //REC::Scintillator
    public void read_Scint_Bank(Bank Scint_Bank){
      if(show){
          System.out.println("REC::Scintillator");
          Scint_Bank.show();
      }
      //Scint_Bank.show();
      for (int k = 0; k < Scint_Bank.getRows(); k++) {
          int pindex = Scint_Bank.getInt("pindex", k);
          double energy = Scint_Bank.getFloat("energy", k);
          int sector = Scint_Bank.getInt("sector", k);
          int layer=Scint_Bank.getInt("layer",k);
          int component=Scint_Bank.getInt("component",k);
          double path = Scint_Bank.getFloat("path", k);
          if (pindex == PIndex && sector==Sector) {
              if(layer==2){
                FTOF_comp=component;
                FTOF_path=path;
                CF_out[9]=component;
                CF_out[10]=(float)path;
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
            if (pindex == PIndex && sector==Sector) {
                Cal_Sector=sector;
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

    public int nCalHits(){
      int nCalHits=0;
      if(PCAL_energy_fcf>0.01){
        nCalHits++;
      }
      if(ECIN_energy_fcf>0.01){
        nCalHits++;
      }
      if(ECOUT_energy_fcf>0.01){
        nCalHits++;
      }

      return nCalHits;
    }

    public void find_RECParticle_fromtrack(Bank RECP,Bank RECT,Bank TTracks, Bank TClust,float[] predtrack,double maxdist){

    
        //System.out.println("\n New pred track");
        //System.out.println(Arrays.toString(predtrack));
        

        for (int i = 0; i < RECP.getRows(); i++) {
            Level3Particle part = new Level3Particle();
            part.read_Particle_Bank(i, RECP);
            if(part.PIndex!=-1){
                part.find_sector_track(RECT);
                part.find_track_clusterIDs(TTracks);
                part.find_track_clusters(TClust);

                
                //System.out.println("true track");
                //System.out.println(Arrays.toString(part.track_clusters));

                float dist=0;
                for (int k=0;k<6;k++){
                    dist+=(part.track_clusters[k]-predtrack[k])*(part.track_clusters[k]-predtrack[k]);
                }

                dist=(float) Math.sqrt(dist);
                //System.out.printf("dist %f max %f\n",dist,maxdist);
                if(dist<maxdist){
                    //System.out.println("Matched!");
                    read_Particle_Bank(i, RECP);
                    PIndex=i;
                    unmatched=false;
                    track_clusters=part.track_clusters;
                    Sector=part.Sector;
                    break;
                }


            }
        }
        
    }

    public void find_ClosestRECParticle_fromPredP(Bank PartBank,double lx, double ly, double lz){
        double min_resPx=9999;
        double min_resPy=9999;
        double min_resPz=9999;
        int best_ind=-1;
        for (int i = 0; i < PartBank.getRows(); i++) {
            read_Particle_Bank(i, PartBank);
            double resPx=Math.abs(Px-Pred_Px);
            double resPy=Math.abs(Py-Pred_Py);
            double resPz=Math.abs(Pz-Pred_Pz);
            if(resPx<min_resPx && resPy<min_resPy && resPz<min_resPz && Pred_Charge==Charge){//
                min_resPx=resPx;
                min_resPy=resPy;
                min_resPz=resPz;
                best_ind=i;
            }
        }
        read_Particle_Bank(best_ind, PartBank);
        if(min_resPx<lx && min_resPy<ly && min_resPz<lz && Vz<12 && Vz>-13 && Math.abs(chi2pid)<5){ //inbending
        //if(min_resPx<lx && min_resPy<ly && min_resPz<lz && Vz<10 && Vz>-18 && Math.abs(chi2pid)<20){ //outbending
            PIndex=best_ind;
            unmatched=false;
        } else{
            PIndex=-1;
            PID=-1;
            unmatched=true;
        }
        
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
        double c2p = PartBank.getFloat("chi2pid", pindex);
        double px = PartBank.getFloat("px", pindex);
        double py = PartBank.getFloat("py", pindex);
        double pz = PartBank.getFloat("pz", pindex);
        double vx = PartBank.getFloat("vx", pindex);
        double vy = PartBank.getFloat("vy", pindex);
        double vz = PartBank.getFloat("vz", pindex);
        if (Math.abs(status) >= 2000 && Math.abs(status) < 4000) {
            PID = pid;
            Charge = charge;
            chi2pid=c2p;
            Px = px;
            Py = py;
            Pz = pz;
            Vx = vx;
            Vy = vy;
            Vz = vz;
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
        
        if(PID==11){
            out_ind=0;
        } else {
            out_ind=1;
        }

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
        csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,", PCALDU_fcf/7.0, PCALDV_fcf/7.0, PCALDW_fcf/7.0)); //loop over at most 7 strips
        csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,", ECINDU_fcf/7.0, ECINDV_fcf/7.0, ECINDW_fcf/7.0)); //3 each side of pred
        csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,", ECOUTDU_fcf/7.0, ECOUTDV_fcf/7.0, ECOUTDW_fcf/7.0));
        
        for (float value : track_clusters) {
            csvLineBuilder.append(String.format("%.6f,", value));
        }
        
        for (double value : HTCC_adcs) {
            csvLineBuilder.append(String.format("%.6f,",value/35000));
        }

        //use just one number for HTCC rather than mostly 0 array
        /*float HTCC_val=0;
        for (double value : HTCC_adcs) {
            HTCC_val+=value;
        }
        HTCC_val=HTCC_val/50000;
        csvLineBuilder.append(String.format("%.6f,",HTCC_val));

        //add Momentum theta phi for pred
        csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,", Pred_P/10, Pred_Theta, (Pred_Phi+3.5)/7.0));*/
        
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
    public String get_csv_out_fromCFPred_noTrack(){
      StringBuilder csvLineBuilder = new StringBuilder();

      csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,", PCAL_energy_fcf/150000.0, ECIN_energy_fcf/150000.0, ECOUT_energy_fcf/150000.0));// /3
      csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,", PCALLU_fcf/500.0, PCALLV_fcf/500.0, PCALLW_fcf/500.0)); // /2000
      csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,", ECINLU_fcf/500.0, ECINLV_fcf/500.0, ECINLW_fcf/500.0));
      csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,", ECOUTLU_fcf/500.0, ECOUTLV_fcf/500.0, ECOUTLW_fcf/500.0));
      csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,", PCALDU_fcf/7.0, PCALDV_fcf/7.0, PCALDW_fcf/7.0)); //loop over at most 7 strips
      csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,", ECINDU_fcf/7.0, ECINDV_fcf/7.0, ECINDW_fcf/7.0)); //3 each side of pred
      csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,", ECOUTDU_fcf/7.0, ECOUTDV_fcf/7.0, ECOUTDW_fcf/7.0));
      
      for (double value : HTCC_adcs) {
        csvLineBuilder.append(String.format("%.6f,",value/35000));
    }

      //add Momentum theta phi for pred
      csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,", Pred_P/10, Pred_Theta, (Pred_Phi+3.5)/7.0));
      
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

    //35 floating-point values and 6 integer values...
    public float[] get_vars_forpid_noTrack(){
      float[] vars_for_pid = new float[32];
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
      
      for (float value : HTCC_adcs) {
          vars_for_pid[n]=value/35000;
          n++;
      }

      //add Momentum theta phi for pred
      vars_for_pid[n]=Pred_P/10;
      n++;

      vars_for_pid[n]=Pred_Theta;
      n++;

      vars_for_pid[n]=(Pred_Phi+(float)3.5)/((float)7.);
      n++;
      
      return vars_for_pid;
  }

    public String get_csv_cf_out(){
        StringBuilder csvLineBuilder = new StringBuilder();
        // Append values from track_clusters array to the StringBuilder
        for (float value : track_clusters) {
            csvLineBuilder.append(String.format("%.6f,", value));
        }

        // 9 first values of CF_out are the calorimeter Ls
        for (int i=0;i<9;i++) {
            csvLineBuilder.append(String.format("%.6f,", CF_out[i]/500)); //500 using LU/LV/LW
        }
        //next two values are FTOF comp and path
        csvLineBuilder.append(String.format("%.6f,", FTOF_comp/62));
        csvLineBuilder.append(String.format("%.6f,", FTOF_path/1000));
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
