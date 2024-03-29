
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
public class Level3Converter {
    NeuralClassifierModel trackfinder = new NeuralClassifierModel();
    EJMLModel cf;
    DataList LtoSconv;

    List<float[]> negTracks = new ArrayList<>();
    List<float[]> posTracks = new ArrayList<>();
    List<float[]> all_permutations = new ArrayList<>();

    public Level3Converter() {

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

    public Boolean arrayNotEmpty(float [] arr){
        for (double value : arr) {
            if(value!=0){
                return true;
            }
        }
        return false;
    }

    public Boolean arrayContains(float [] arr,int c){
        for (double value : arr) {
            if(value==c){
                return true;
            }
        }
        return false;
    }

    public void find_tracks_fromClusters(Bank ClusterBank, List<float[]> negTracks, List<float[]> posTracks) {

        //ClusterBank.show();
        all_permutations.clear();

        for (int sect=1;sect<7;sect++){
            List<Float> SL1 = new ArrayList<>();
            List<Float> SL2 = new ArrayList<>();
            List<Float> SL3 = new ArrayList<>();
            List<Float> SL4 = new ArrayList<>();
            List<Float> SL5 = new ArrayList<>();
            List<Float> SL6 = new ArrayList<>();

            for(int i=0;i<ClusterBank.getRows();i++){
                int sl=ClusterBank.getInt("superlayer", i);
                int s=ClusterBank.getInt("sector",i);
                float avgw=ClusterBank.getFloat("avgWire", i)/112;
                //System.out.println("Adding sls");
                if(s==sect){
                    if(sl==1){
                        SL1.add(avgw);
                    } else if(sl==2){
                        SL2.add(avgw);
                    } else if(sl==3){
                        SL3.add(avgw);
                    }else if(sl==4){
                        SL4.add(avgw);
                    }else if(sl==5){
                        SL5.add(avgw);
                    }else if(sl==6){
                        SL6.add(avgw);
                    }
                }
            }

            // Combine all SLs into a List of Lists
            List<List<Float>> SLs = Arrays.asList(
                    SL1, SL2, SL3, SL4, SL5, SL6);

            // Create a list to store permutations
            List<float[]> permutations = new ArrayList<>();

            //System.out.println("gen perms");

            // Call the generatePermutations method to populate the permutations list
            generatePermutations(SLs, 0, new float[6], permutations,50.0f / 112.0f);

            //for debugging
            for(float[] p:permutations){all_permutations.add(p);}


            //System.out.println("Adding perms");
            // Print the permutations
            for (float[] permutation : permutations) {
                float[] tr_pred= new float[3];
                trackfinder.getModel().getOutput(permutation,tr_pred);
                if(tr_pred[2]>0.1){
                    negTracks.add(permutation);
                } else if(tr_pred[1]>0.1){
                    posTracks.add(permutation);
                }

            }

        }
        
        
    }

    // Recursive method to generate permutations with constraint
    private static void generatePermutations(List<List<Float>> vectors, int depth, float[] current, List<float[]> permutations, float maxDifference) {
        if (depth == vectors.size()) {
            // Base case: we've filled in all slots, add this permutation to the list
            permutations.add(Arrays.copyOf(current, current.length));
            return;
        }

        // Iterate over the elements in the current vector
        for (Float element : vectors.get(depth)) {
            // Check the constraint with the previous element
            if (depth > 0 && Math.abs(element - current[depth - 1]) > maxDifference) {
                continue;  // Skip this element if the constraint is violated
            }

            // Place the current element in the current slot
            current[depth] = element;

            // Recursively generate permutations for the next depth
            generatePermutations(vectors, depth + 1, current, permutations, maxDifference);
        }
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

    public int nPart_pSect(List<Level3Candidate> Candidates,int sect){
        int nPart_pSect=0;
        for(Level3Candidate part:Candidates){
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

    //                              0           1           2               3           4
    //                  5               6                   7               8               9            
    //                  10                      11
    //Bank[] banks = r.getBanks("DC::tdc", "ECAL::adc", "RUN::config", "FTOF::adc", "HTCC::adc",
    //            "REC::Particle", "REC::Calorimeter", "REC::Cherenkov","REC::Track","HitBasedTrkg::HBClusters",
    //            "HitBasedTrkg::HBTracks","ECAL::clusters");
    public List<Level3Candidate> getCandidates_ownc(Bank[] banks){
        

        negTracks.clear();
        posTracks.clear();
        
        find_tracks_fromClusters(banks[9],negTracks,posTracks);

        //Instant start_time = Instant.now();

        List<Level3Candidate> ps = new ArrayList<Level3Candidate>();
        for(float[] negtrack:negTracks){
        // find and initialise Candidates
            Level3Candidate part = new Level3Candidate();
            part.setShow(true);

            part.find_RECParticle_fromtrack(banks[5],banks[8],banks[10],banks[9],negtrack,0.05);
            part.read_Cal_Bank(banks[6]);
            part.find_HTCC_ADCs(banks[4]);
            part.set_pid_label();
            //part.print();
            //part.printbanks(banks[5],banks[6],banks[4]);
            ps.add(part);
        }


        //Instant finish_time = Instant.now();
        //timeElapsed += Duration.between(start_time, finish_time).toMillis();

        return ps;
    }

    //                              0           1           2               3           4
    //                  5               6                   7               8               9            
    //                  10                      11
    //Bank[] banks = r.getBanks("DC::tdc", "ECAL::adc", "RUN::config", "FTOF::adc", "HTCC::adc",
    //            "REC::Particle", "REC::Calorimeter", "REC::Cherenkov","REC::Track","HitBasedTrkg::HBClusters",
    //            "HitBasedTrkg::HBTracks","ECAL::clusters");
    public List<Level3Candidate> getCandidates(Bank[] banks,CompositeNode tr,CompositeNode pt){

        //Instant start_time = Instant.now();

        List<Level3Candidate> ps = new ArrayList<Level3Candidate>();
        for(int i = 0; i < tr.getRows(); i++){
        // find and initialise Candidates
            Level3Candidate part = new Level3Candidate();
            //part.setShow(true);
            part.find_AI_tracks(tr,pt,i);
            part.find_ClosestRECParticle_fromPredP(banks[5],999,999,999);
            part.read_Cal_Bank(banks[6]);
            part.find_HTCC_ADCs(banks[4]);
            part.set_pid_label();

            float[] cf_pred= new float[9];
            cf.getOutput(part.track_clusters, cf_pred);
            int[] cf_strips=convLtoStrip(cf_pred);
            part.read_cal_bank_from_cf_pred(cf_pred,cf_strips,banks[1]);

            //part.print();
            //part.printbanks(tr,pt,banks[5],banks[6],banks[4]);

            ps.add(part);
        }

        //Instant finish_time = Instant.now();
        //timeElapsed += Duration.between(start_time, finish_time).toMillis();

        return ps;
    }

    public int hasMatchedPID(List<Level3Candidate> Candidates,int pid){
        int c = 0;
        for (Level3Candidate part : Candidates) {
            int PID=part.PID;
            if (PID == pid) {
                c++;
            }
        }
        return c;
    }

    public int eventHasRECPID(Bank[] banks,int PID,ArrayList<Level3Particle> parts,CompositeNode tr,CompositeNode pt){
        int c = 0;
        for (int i = 0; i < banks[5].getRows(); i++) {
            Level3Particle part = new Level3Particle();
            part.read_Particle_Bank(i, banks[5]);
            part.find_sector_track(banks[8]);
            part.find_track_clusterIDs(banks[10]);
            part.find_track_clusters(banks[9]);
            part.read_Cal_Bank(banks[6]);
            part.find_sector_cal(banks[6]);
            part.find_HTCC_ADCs(banks[4]);
            
            if (part.Sector!=0 && part.Track_nSL==6 && part.Track_chi2<350) {
                parts.add(part);

                //part.applyTriangCut();
                if(PID == part.PID){
                    if(PID==11 ){ //&& arrayNotEmpty(part.HTCC_adcs)
                        if(part.Vz<12 && part.Vz>-13 && Math.abs(part.chi2pid)<5){
                            c++;
                        }
                    } else{
                        c++;
                    }
                    /*tr.print();
                    pt.print();
                    banks[5].show();*/
                    
                }
                
            }
        }
        return c;
    }

    public void convertData(String file,String out,Boolean isPID,int nPart){

        // Specify the path where you want to save the CSV file
        Path filePath = Paths.get(out);

        // Delete the file if it already exists
        try {
            Files.deleteIfExists(filePath);
        } catch (IOException e) {
            e.printStackTrace();
        }

        Boolean notAllFull=true;
        int nEl=0, nPos=0, nPip=0, nPim=0,nMum=0,nMup=0,rEv=0;
        int nUnnmatchedEl=0;

        try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(filePath))) {
            // String file=dir;

            HipoReader r = new HipoReader(file);

            Event e = new Event();

            // r.getSchemaFactory().show();
            Bank[] banks = r.getBanks("DC::tdc", "ECAL::adc", "RUN::config", "FTOF::adc", "HTCC::adc",
                "REC::Particle", "REC::Calorimeter", "REC::Cherenkov","REC::Track","HitBasedTrkg::HBClusters",
                "HitBasedTrkg::HBTracks","ECAL::clusters");

            CompositeNode tr = new CompositeNode(3210,2,"i",1200);
            CompositeNode pt = new CompositeNode(3210,3,"i",1200);

            while (r.hasNext() && notAllFull) {
                rEv++;

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

                e.read(tr,32100,2);
                e.read(pt,32100,3);
    
                List<Level3Candidate> Candidates = getCandidates(banks,tr,pt);
                //List<Level3Candidate> Candidates = getCandidates_ownc(banks);
                for (Level3Candidate p : Candidates) {
                    if (p.unmatched == false && p.Cal_Sector == p.Sector && p.Charge == -1) {
                        //p.applyTriangCut();
                        if (p.PID == 11 && nEl < nPart) {
                            if (isPID) {
                                if (p.PCALLU_fcf != -2) { // && arrayNotEmpty(p.HTCC_adcs)
                                    writer.println(p.get_csv_out_fromCFPred());
                                    nEl++;
                                }
                                // writer.println(p.get_csv_out());
                                // nEl++;
                            } else {
                                if (!arrayContains(p.CF_out, 0)) {
                                    writer.println(p.get_csv_cf_out());
                                    nEl++;
                                }
                            }
                        } /*
                           * else if (p.PID == -211 && nPim < nPart ) {
                           * if (isPID) {
                           * if (p.PCALLU_fcf != -2) {
                           * writer.println(p.get_csv_out_fromCFPred());
                           * nPim++;
                           * }
                           * // writer.println(p.get_csv_out());
                           * // nPim++;
                           * }
                           * }
                           */
                        // any neg particle not an e- as BG
                        else if (p.PID != 11 && nPim < nPart) {
                            if (isPID) {
                                if (p.PCALLU_fcf != -2) {
                                    writer.println(p.get_csv_out_fromCFPred());
                                    nPim++;
                                }
                                // writer.println(p.get_csv_out());
                                // nPim++;
                            }
                        }
                    }
                }

                ArrayList<Level3Particle> parts=new ArrayList<Level3Particle>();
                int nRecEl=eventHasRECPID(banks, 11,parts,tr,pt);
                int nMatchedEl=hasMatchedPID(Candidates, 11);
                if(nRecEl>nMatchedEl){
                    /*System.out.printf("\nnREC e- %d nMatched e- %d \n",nRecEl,nMatchedEl);
                    System.out.println("AI tracks");
                    for(int i=0;i<tr.getRows();i++){
                        System.out.printf("sector = %2d, charge = %3d segments [ ", tr.getInt(1,i),tr.getInt(2,i));
                        for(int k = 0; k < 6; k++) System.out.printf(" %9.5f ",tr.getDouble(10+k,i));
                        System.out.printf("] p %9.5f  %9.5f %9.5f\n",pt.getDouble(4,i),pt.getDouble(5,i),pt.getDouble(6,i)); 
                    }
                    System.out.println("REC tracks");
                    for (Level3Particle p:parts){
                        if(p.Charge!=0){
                            System.out.printf("%d Sector %d charge %d segments [ ",p.PID,p.Sector, p.Charge);
                            for(int k = 0; k < 6; k++) System.out.printf(" %9.5f ",p.track_clusters[k]*112);
                            System.out.printf("] p %f %f %f\n",p.Px,p.Py,p.Pz);
                        }
                        
                    }
                    banks[5].show();
                    banks[8].show();*/
                    //banks[9].show();

                    /*System.out.println("\nall possible tracks");
                    for (float[] ntr:all_permutations){
                        System.out.println(Arrays.toString(ntr));
                    }
                    System.out.println("\nneg pred tracks");
                    for (float[] ntr:negTracks){
                        System.out.println(Arrays.toString(ntr));
                    }
                    System.out.println("pos pred tracks");
                    for (float[] ptr:posTracks){
                        System.out.println(Arrays.toString(ptr));
                    }
                    System.out.println("parts");
                    for (Level3Particle p:parts){
                        System.out.printf("part PID %d Pindex %d Sector %d", p.PID, p.PIndex, p.Sector);
                        System.out.println(Arrays.toString(p.track_clusters));
                    }

                    System.out.println("REC::Track");
                    banks[8].show();
                    //System.out.println("HB::Cluster");
                    //banks[9].show();
                    System.out.println("HB::Tracks");
                    banks[10].show();*/

                    nUnnmatchedEl+=nRecEl-nMatchedEl;
                }

                if (nEl >= nPart && nPim >= nPart) { // && nPip>=nPart ){// && nPos>= nPart && nMum>=nPart && nMup>=
                                                     // nPart){
                    notAllFull = false;
                }

            }

        } catch (IOException e) {
            e.printStackTrace();
        }

        double eff=((float) nEl)/(((float) nEl)+((float) nUnnmatchedEl));
        System.out.printf("Found %d 11, %d -11, %d -211, %d 211, %d -13, %d 13\n",nEl,nPos,nPim,nPip,nMum,nMup);
        System.out.printf("Nb UnMatched 11: %d, eff %f\n",nUnnmatchedEl,eff);
        System.out.printf("Read %d events\n\n",rEv);

    }

    
    
    public static void main(String[] args){        

        //to run, in order
        // /open Level3Particle.java
        // /open Level3Candidate.java
        // /open Level3Converter.java
        // Level3Converter.main(new String[]{});

        //String file = "/Users/tyson/data_repo/trigger_data/rgd/018326/recook_caos_pid/run_018326_3_wAIBanks.h5";
        //String filet = "/Users/tyson/data_repo/trigger_data/rgd/018326/recook_caos_pid/run_018326_4_wAIBanks.h5";
        //String outDir="/Users/tyson/data_repo/trigger_data/rgd/018326/for_caos_pid/";

        /*String file = "/Users/tyson/data_repo/trigger_data/rgd/018640/recook_wAI/run_018640_1_wAI.h5";
        String filet = "/Users/tyson/data_repo/trigger_data/rgd/018640/recook_wAI/run_018640_2_wAI.h5";
        String outDir = "/Users/tyson/data_repo/trigger_data/rgd/018640/for_caos_pid/";*/

        String file = "/Users/tyson/data_repo/trigger_data/rgd/018437_AI/rec_clas_018437.evio.00000-00004.hipo";
        String filet = "/Users/tyson/data_repo/trigger_data/rgd/018437_AI/rec_clas_018437.evio.00005-00009.hipo";
        String outDir = "/Users/tyson/data_repo/trigger_data/rgd/018437_AI/for_caos_pid/";

        Level3Converter conv = new Level3Converter();
        conv.load_trackfinder("clas12rgd.network");
        conv.load_cf("cf_el.network");
        conv.load_LtoStripConv("LtoStrip_convTable.csv");
        
        //conv.convertData(file,outDir+"train_fromcfpred_allNegBG2.csv",true,75000);//50000
        //conv.convertData(filet,outDir+"test_fromcfpred_allNegBG2.csv",true,10000);//50000

        conv.convertData(file,outDir+"train_cf.csv",false,100000);//50000
        conv.convertData(filet,outDir+"test_cf.csv",false,50000);//50000

        //conv.convertData(filet,outDir+"bla_cf.csv",false,91,90,100);//50000
        
        
       

    }
}
