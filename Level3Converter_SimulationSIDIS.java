
import j4np.data.base.DataFrame;
import j4np.hipo5.data.Bank;
import j4np.hipo5.data.CompositeNode;
import j4np.hipo5.data.Event;
import j4np.hipo5.data.Node;
import j4np.hipo5.io.HipoReader;
import j4np.hipo5.io.HipoWriter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import java.util.List;
import java.io.File;
import java.io.IOException;
import java.lang.Math;

import twig.data.GraphErrors;
import twig.data.H1F;
import twig.data.H2F;
import twig.graphics.TGCanvas;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.io.IOException;
import java.io.PrintWriter;

import j4ml.data.*;
import j4ml.deepnetts.*;
import j4ml.ejml.EJMLModel;
import j4ml.ejml.EJMLModel.ModelType;
import j4np.neural.classifier.NeuralClassifierModel;

/**
 *
 * @author tyson
 */
public class Level3Converter_SimulationSIDIS {
    NeuralClassifierModel trackfinder = new NeuralClassifierModel();
    EJMLModel cf;
    DataList LtoSconv;

    public Level3Converter_SimulationSIDIS() {

    }

    public void load_trackfinder(String path){
        trackfinder.loadFromFile(path,12);
    }

    public void load_cf(String path){
        cf = new EJMLModel(path, ModelType.TANH_LINEAR);
    }

    public void load_LtoStripConv(String path){
        LtoSconv = DataList.fromCSV("LtoStrip_convTable.csv",
                DataList.range(0,69), DataList.range(0,1));
    }

    public Boolean arrayContains(float [] arr,int c){
        for (double value : arr) {
            if(value==c){
                return true;
            }
        }
        return false;
    }

    int[] convLtoStrip(float[] Ls){
        int[] strips = new int[]{-2,-2,-2,-2,-2,-2,-2,-2,-2};
        int[] nStrips = new int[]{68,63,63,37,37,37,37,37,37};
        for(int i=0;i<9;i++){
            Ls[i]=Ls[i]*500;
            float[] conv_det_view = LtoSconv.getList().get(i).floatFirst();
            for (int str=0;str<(nStrips[i]);str++){
                //in PCAL V, W take distance from other end
                //System.out.printf("wtf %d %d \n",i,str);
                if(i==1 || i==2){
                    if(Ls[i]<conv_det_view[str] && Ls[i]>=conv_det_view[str+1]){
                        strips[i]=str+1; //add 1 as the upper limit applies to strip+1 due to inverted order
                    }
                } else {
                    if(Ls[i]>=conv_det_view[str] && Ls[i]<conv_det_view[str+1]){
                        strips[i]=str;
                    }
                }
            }
        }
        return strips;
    }

    public int nPart_pSect(List<Level3Particle> particles,int sect){
        int nPart_pSect=0;
        for(Level3Particle part:particles){
            if(part.Cal_Sector==sect){
                nPart_pSect++;
            }
        }
        return nPart_pSect;
    }

    public int nTrack_pSect(Bank TrackBank,int sect){
        int nTrack_pSect=0;
        for (int k = 0; k < TrackBank.getRows(); k++) {
            int pindex = TrackBank.getInt("pindex", k);
            int sectorTrk = TrackBank.getInt("sector", k);
            if(sectorTrk==sect){nTrack_pSect++;}
        }
        return nTrack_pSect;
    }

    /*                          //0         1           2               3           4
    //  5               6               7                   8
    //  9               10              11                      12
    //  13                              14                      15
    Bank[] banks = r.getBanks("DC::tdc", "ECAL::adc", "RUN::config", "FTOF::adc", "HTCC::adc",
    "REC::Particle", "REC::Track", "REC::Calorimeter", "REC::Cherenkov",
    "ECAL::clusters", "MC::Particle","HitBasedTrkg::HBTracks","HitBasedTrkg::HBClusters",
    "TimeBasedTrkg::TBTracks","TimeBasedTrkg::TBClusters","ECAL::hits");*/

    public List<Level3Particle> getRECPart(Bank[] dsts,int nEl){
        List<Level3Particle> ps = new ArrayList<Level3Particle>();

        // find and initialise particles
        for (int i = 0; i < dsts[5].getRows(); i++) {
            Level3Particle part = new Level3Particle();
            part.read_Particle_Bank(i, dsts[5]);
            if(part.PIndex!=-1){ //ie part in FD
                part.find_ClosestMCParticle(dsts[10]);
                //for debugging
                /*if (part.MC_PID == 11 && part.TruthMatch(0.1, 0.1, 0.1)&&nEl==0) {
                    part.setShow(true);
                    System.out.printf("\n New Particle pindex %d \n",i);
                    System.out.println("REC::Particle");
                    dsts[5].show();
                }*/
                part.read_Cal_Bank(dsts[7]);
                part.read_HTCC_bank(dsts[8]);
                part.find_sector_cal(dsts[7]);
                part.find_sector_track(dsts[6]);
                part.find_track_clusterIDs(dsts[11]);
                part.find_track_clusters(dsts[12]);
                part.find_HTCC_ADCs(dsts[4]);
                part.set_pid_label();

                float[] tr_pred= new float[3];
                trackfinder.getModel().getOutput(part.track_clusters,tr_pred);
                
                //for this model pred_1 is negatives
                if(tr_pred[1]>0.9){
                    float[] cf_pred= new float[9];
                    cf.getOutput(part.track_clusters, cf_pred);
                    int[] cf_strips=convLtoStrip(cf_pred);
                    part.read_cal_bank_from_cf_pred(cf_pred,cf_strips,dsts[1]);
                }

                ps.add(part);
                
            }
        }
        return ps;
    }

    public int hasPID(List<Level3Particle> particles,int pid,boolean REC){
        int c = 0,el_i=-1;
        for (Level3Particle part : particles) {
            int PID=part.PID;
            if(!REC){
                PID=part.MC_PID;
            } 
            if (PID == pid) {
                el_i = c;
            }
            c++;
        }
        return el_i;
    }

    public void convertData(String dir,String out,Boolean isPID,int nFiles,int start,int nPart){

        // Specify the path where you want to save the CSV file
        Path filePath = Paths.get(out);

        // Delete the file if it already exists
        try {
            Files.deleteIfExists(filePath);
        } catch (IOException e) {
            e.printStackTrace();
        }

        Boolean notAllFull=true;
        int nEl=0, nPos=0, nPip=0, nPim=0,nMum=0,nMup=0;

        try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(filePath))) {
            while (start < nFiles && notAllFull) {// 56

                String file = dir + "clasdis_" + String.valueOf(start) + ".hipo";
                // String file=dir;
                start++;

                HipoReader r = new HipoReader(file);

                Event e = new Event();

                // r.getSchemaFactory().show();
                Bank[] banks = r.getBanks("DC::tdc", "ECAL::adc", "RUN::config", "FTOF::adc", "HTCC::adc",
                        "REC::Particle", "REC::Track", "REC::Calorimeter", "REC::Cherenkov",
                        "ECAL::clusters", "MC::Particle", "HitBasedTrkg::HBTracks", "HitBasedTrkg::HBClusters",
                        "TimeBasedTrkg::TBTracks", "TimeBasedTrkg::TBClusters", "ECAL::hits");

                while (r.hasNext() && notAllFull) {

                    r.nextEvent(e);
                    e.read(banks);

                    // System.out.println("DC");
                    // banks[0].show();

                    /*
                     * System.out.println("FTOF");
                     * banks[3].show();
                     * System.out.println("HTCC");
                     * banks[4].show();
                     */

                     /*
                     * System.out.println("REC::Calorimeter");
                     * banks[7].show();
                     * System.out.println("ECAL::hits");
                     * banks[15].show();
                     */

                    List<Level3Particle> particles = getRECPart(banks,nEl);
                    for (Level3Particle p : particles) {
                        if (p.TruthMatch(0.1, 0.1, 0.1) && p.ECAL_energy!=0) {
                            if (p.MC_PID == 11 && nEl < nPart && p.HTCC_Sector==p.Sector && p.Cal_Sector==p.Sector) {
                                if(isPID){
                                    if(p.PCALLU_fcf!=-2){
                                        writer.println(p.get_csv_out_fromCFPred());
                                        nEl++;
                                    }
                                    //writer.println(p.get_csv_out());
                                    //nEl++;
                                } else{
                                    if(p.ECIN_energy!=0 && !arrayContains(p.CF_out,0)){
                                        writer.println(p.get_csv_cf_out());
                                        nEl++;
                                    }
                                    
                                }
                                
                            } else if (p.MC_PID == -211 && nPim < nPart && p.Cal_Sector==p.Sector) {
                                if(isPID){
                                    if(p.PCALLU_fcf!=-2){
                                        writer.println(p.get_csv_out_fromCFPred());
                                        nPim++;
                                    }
                                    //writer.println(p.get_csv_out());
                                    //nPim++;
                                } /*else{ //only use electrons for cf
                                    if(p.ECIN_energy!=0){
                                        writer.println(p.get_csv_cf_out());
                                        nPim++;
                                    }
                                }*/
                               
                            } /*else if (p.MC_PID == 211 && nPip < nPart && p.Cal_Sector==p.Sector) {
                                if(isPID){
                                    writer.println(p.get_csv_out());
                                } else{
                                    if(p.ECIN_energy!=0){
                                        writer.println(p.get_csv_cf_out());
                                        nPip++;
                                    }
                                }
                                
                            } else if (p.MC_PID == -11 && nPos < nPart && p.HTCC_Sector==p.Sector && p.Cal_Sector==p.Sector) {
                                if(isPID){
                                    writer.println(p.get_csv_out());
                                    nPos++;
                                } else{
                                    if(p.ECIN_energy!=0){
                                        writer.println(p.get_csv_cf_out());
                                        nPos++;
                                    }
                                }
                                
                            } else if (p.MC_PID == -13 && nMum < nPart && p.Cal_Sector==p.Sector) {
                                if(isPID){
                                    writer.println(p.get_csv_out());
                                    nMum++;
                                } else{
                                    if(p.ECIN_energy!=0){
                                        writer.println(p.get_csv_cf_out());
                                       nMum++;
                                    }
                                }
                                
                            } else if (p.MC_PID == 13 && nMup < nPart && p.Cal_Sector==p.Sector) {
                                if(isPID){
                                    writer.println(p.get_csv_out());
                                    nMup++;
                                } else{
                                    if(p.ECIN_energy!=0){
                                        writer.println(p.get_csv_cf_out());
                                       nMup++;
                                    }
                                }
                                
                            }*/
                        }
                    }

                    if(nEl >= nPart && nPim>=nPart){ //&& nPip>=nPart ){// && nPos>= nPart && nMum>=nPart && nMup>= nPart){
                        notAllFull=false;
                    }

                }
                
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.printf("Found %d 11, %d -11, %d -211, %d 211, %d -13, %d 13\n",nEl,nPos,nPim,nPip,nMum,nMup);

    }

    
    
    public static void main(String[] args){        

        //to run
        // /open Level3Particle.java
        // /open Level3Converter_SimulationSIDIS.java
        // Level3Converter_SimulationSIDIS.main(new String[]{});

        String dir = "/Users/tyson/data_repo/trigger_data/sims/claspyth_train/";

        Level3Converter_SimulationSIDIS conv = new Level3Converter_SimulationSIDIS();
        conv.load_trackfinder("clas12rgd.network");
        conv.load_cf("cf_el.network");
        conv.load_LtoStripConv("LtoStrip_convTable.csv");
        
        conv.convertData(dir,dir+"for_pid/train_fromcfpred.csv",true,90,1,98000);//50000
        conv.convertData(dir,dir+"for_pid/test_fromcfpred.csv",true,98,90,8000);//50000
        //conv.convertData(dir,dir+"for_pid/train_cf.csv",false,90,1,1000000);//50000
        //conv.convertData(dir,dir+"for_pid/test_cf.csv",false,98,90,50000);//50000
        //conv.convertData(dir,dir+"for_pid/bla_cf.csv",false,91,90,100);//50000
        
        
       

    }
}
