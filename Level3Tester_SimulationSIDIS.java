
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
import java.time.Duration;
import java.time.Instant;
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
public class Level3Tester_SimulationSIDIS {
    NeuralClassifierModel trackfinder = new NeuralClassifierModel();
    EJMLModel cf;
    EJMLModel pider;
    DataList LtoSconv;
    float timeElapsed=0;

    public Level3Tester_SimulationSIDIS() {

    }

    public void load_trackfinder(String path){
        trackfinder.loadFromFile(path,12);
    }

    public void load_pider(String path){
        pider= new EJMLModel(path, ModelType.SOFTMAX);
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

    public List<Level3Particle> getRECPart(Bank[] dsts){
        //Instant start_time = Instant.now();

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
                part.setTrackResp(tr_pred[1]);
                //for this model pred_1 is negatives
                //if(tr_pred[1]>0.9){
                float[] cf_pred= new float[9];
                cf.getOutput(part.track_clusters, cf_pred);
                int[] cf_strips=convLtoStrip(cf_pred);
                part.read_cal_bank_from_cf_pred(cf_pred,cf_strips,dsts[1]);
                float[] pid_pred =new float[2];
                pider.getOutput(part.get_vars_forpid(), pid_pred);
                part.setPidResp(pid_pred[0]);
                    
                //}

                /*if (part.MC_PID == 11 && part.TruthMatch(0.1, 0.1, 0.1)) {
                    System.out.println(Arrays.toString(part.get_vars_forpid()));
                    System.out.printf("tr_pred %f cf pred %d pid_pred %f \n", tr_pred[1],cf_strips[0],pid_pred[0]);
                }*/

                ps.add(part);
                
            }
        }


        //Instant finish_time = Instant.now();
        //timeElapsed += Duration.between(start_time, finish_time).toMillis();

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

    public DataList getPred(String dir,int nFiles,int start, int nPred){

        

        DataList pred = new DataList();

        Boolean notAllFull=true;
        int nRead=0;
        float nEv=0;

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

            Instant start_time = Instant.now();
            while (r.hasNext() && notAllFull) {
                nEv++;

                r.nextEvent(e);
                e.read(banks);

                List<Level3Particle> particles = getRECPart(banks);
                for (Level3Particle p : particles) {
                    if (p.TruthMatch(0.1, 0.1, 0.1) && p.track_resp > 0.9) {
                        double phi=p.Phi*(180.0 / Math.PI);
                        double theta=p.Theta*(180.0 / Math.PI);
                        int mcpid=p.MC_PID;
                        if(mcpid==11){
                            /*if (p.HTCC_Sector != p.Sector) {
                                mcpid = 0;
                            }*/
                        }
                        DataEntry de = new DataEntry(new double[]{mcpid,p.P,theta,phi,p.track_resp}, new double[]{p.pid_resp});
                        pred.add(de);
                        nRead++;
                    }
                }

                if (nRead>=nPred) { 
                    notAllFull = false;
                }
            } 
            Instant finish_time = Instant.now();
            System.out.println("t s"+start_time.toString()+" f "+finish_time.toString());
            timeElapsed += Duration.between(start_time, finish_time).toMillis();
        }
        System.out.printf("Found %d particles\n",nRead);

        float rate=(nEv/timeElapsed);
        System.out.printf("Took %f ms for %f events ie rate %f kHz\n\n",timeElapsed,nEv,rate);

        return pred;

    }

    public void PlotResponse(DataList pred, int elClass, int desiredPID,String part) {
        int NEvents = pred.getList().size();
    
        H1F hRespPos = new H1F(part+" in Sector", 100, 0, 1);
        hRespPos.attr().setLineColor(2);
        hRespPos.attr().setFillColor(2);
        hRespPos.attr().setLineWidth(3);
        hRespPos.attr().setTitleX("Response");
        H1F hRespNeg = new H1F("No "+part+" in Sector", 100, 0, 1);
        hRespNeg.attr().setLineColor(5);
        hRespNeg.attr().setLineWidth(3);
        hRespNeg.attr().setTitleX("Response");
        //Sort predictions into those made on the positive/or negative samples
        for(int i=0;i<NEvents;i+=1) {
            float[] vars = pred.getList().get(i).floatFirst();
            float[] resp = pred.getList().get(i).floatSecond();
    
            if(vars[0]==desiredPID) {
                hRespPos.fill(resp[elClass]);
            } else {
                hRespNeg.fill(resp[elClass]);
            }
        }
    
        TGCanvas c = new TGCanvas();
        
        c.setTitle("Response");
        c.draw(hRespPos).draw(hRespNeg,"same");
        c.region().showLegend(0.05, 0.95);
            
    }//End of PlotResponse

    public void PlotVar(DataList pred, int varIndex, int desiredPID,String part,String varName, String varUnits,double low, double high,int nBins) {
        int NEvents = pred.getList().size();
    
        H1F hRespPos = new H1F(part+" in Sector", nBins,low,high);
        hRespPos.attr().setLineColor(2);
        hRespPos.attr().setFillColor(2);
        hRespPos.attr().setLineWidth(3);
        hRespPos.attr().setTitleX(varName+" "+varUnits);
        H1F hRespNeg = new H1F("No "+part+" in Sector", nBins,low,high);
        hRespNeg.attr().setLineColor(5);
        hRespNeg.attr().setLineWidth(3);
        hRespNeg.attr().setTitleX(varName+" "+varUnits);
        //Sort predictions into those made on the positive/or negative samples
        for(int i=0;i<NEvents;i+=1) {
            float[] vars = pred.getList().get(i).floatFirst();
    
            if(vars[0]==desiredPID) {
                hRespPos.fill(vars[varIndex]);
            } else {
                hRespNeg.fill(vars[varIndex]);
            }
        }
    
        TGCanvas c = new TGCanvas();
        
        c.setTitle(varName);
        c.draw(hRespPos).draw(hRespNeg,"same");
        c.region().showLegend(0.05, 0.95);
            
    }//End of PlotResponse
    
    
    
    public double[] getMetrics(DataList pred, int elClass, int desiredPID,double thresh){
        double[] metrics= new double[5];
        int nEvents = pred.getList().size();
    
        double TP=0,FP=0,FN=0;
        for (int i = 0; i < nEvents; i++) {
            float[] vars = pred.getList().get(i).floatFirst();
            float[] resp = pred.getList().get(i).floatSecond();
            if (vars[0]==desiredPID) {
                if (resp[elClass] > thresh) {
                    TP++;
                } else {
                    FN++;
                } 
            } else {
                if (resp[elClass] > thresh) {
                    FP++;
                } 
            } // Check true label
        }
        double Pur=TP/(TP+FP);
        double Eff=TP/(TP+FN);
        metrics[0]=Pur;
        metrics[1]=Eff;
        metrics[2]=TP;
        metrics[3]=FP;
        metrics[4]=FN;
    
        /*System.out.printf("Theres %d electrons in sample\n", nEls);
        System.out.printf("L1 trigger fired %d times in sample\n", nTrig);*/
        return metrics;
    }
    
    //Labels col 0 is 1 if there's an e-, 0 otherwise
    public double[] getMetsForBin(DataList pred, int elClass, int desiredPID,double thresh,int cutVar,double low,double high){
        double[] metrics = new double [2];
        int nEvents = pred.getList().size();
    
        double TP=0,FN=0,FP=0;
        for (int i = 0; i < nEvents; i++) {
            float[] vars = pred.getList().get(i).floatFirst();
            float[] resp = pred.getList().get(i).floatSecond();
            if (vars[cutVar] > low && vars[cutVar]<high) {
                if (vars[0]==desiredPID) {
                    if (resp[elClass] > thresh) {
                        TP++;
                    } else {
                        FN++;
                    } 
                } else {
                    if (resp[elClass] > thresh) {
                        FP++;
                    } 
                } // Check true label
            }
        }
        double Pur=TP/(TP+FP);
        double Eff=TP/(TP+FN);
        metrics[0]= Pur;
        metrics[1]= Eff;
    
        return metrics;
    }
    
    public void plotVarDep(DataList pred, int elClass, int desiredPID,double thresh,Boolean addPur, 
            int cutVar, String varName, String varUnits,double low, double high,double step) {
    
            String yTitle="Metrics";
            if(!addPur){yTitle="Efficiency";}
    
            GraphErrors gEff = new GraphErrors();
            gEff.attr().setMarkerColor(2);
            gEff.attr().setMarkerSize(10);
            gEff.attr().setTitle("Level3 Efficiency");
            gEff.attr().setTitleX(varName+" "+varUnits);
            gEff.attr().setTitleY(yTitle);
    
            GraphErrors gPur = new GraphErrors();
            gPur.attr().setMarkerColor(5);
            gPur.attr().setMarkerSize(10);
            gPur.attr().setTitle("Level3 Purity");
            gPur.attr().setTitleX(varName+" "+varUnits);
            gPur.attr().setTitleY(yTitle);
    
    
            for (double q2=low;q2<high;q2+=step){
                double[] metrics=getMetsForBin(pred,elClass,desiredPID,thresh,cutVar,q2,q2+step);
                gPur.addPoint(q2+step/2, metrics[0], 0, 0);
                gEff.addPoint(q2+step/2, metrics[1], 0, 0);
            } // Increment threshold on response
    
            
    
            TGCanvas c = new TGCanvas();
            c.setTitle("Efficiency vs "+varName);
            c.draw(gEff);
            if(addPur){c.draw(gPur, "same");}
            c.region().axisLimitsY(gPur.getVectorY().getMin()-0.1, 1.05);
            c.region().showLegend(0.6, 0.25);
            
    
    }
    
    public double findBestThreshold(DataList pred, int elClass, int desiredPID,double effLow){
        
        GraphErrors gEff = new GraphErrors();
        gEff.attr().setMarkerColor(2);
        gEff.attr().setMarkerSize(10);
        gEff.attr().setTitle("Efficiency");
        gEff.attr().setTitleX("Response");
        gEff.attr().setTitleY("Metrics");
        GraphErrors gPur = new GraphErrors();
        gPur.attr().setMarkerColor(5);
        gPur.attr().setMarkerSize(10);
        gPur.attr().setTitle("Purity");
        gPur.attr().setTitleX("Response");
        gPur.attr().setTitleY("Metrics");
        double bestRespTh = 0;
        double bestPuratEffLow= 0;
    
        // Loop over threshold on the response
        for (double RespTh = 0.01; RespTh < 0.99; RespTh += 0.01) {
            double metrics[]=getMetrics(pred,elClass,desiredPID,RespTh);
            double Pur = metrics[0];
            double Eff = metrics[1];
            gPur.addPoint(RespTh, Pur, 0, 0);
            gEff.addPoint(RespTh, Eff, 0, 0);
            if (Eff > effLow) {
                if (Pur > bestPuratEffLow) {
                    bestPuratEffLow = Pur;
                    bestRespTh = RespTh;
                }
            }
        } // Increment threshold on response
    
        System.out.format("%n Best Purity at Efficiency above %f: %.3f at a threshold on the response of %.3f %n%n",
                effLow*100,bestPuratEffLow*100, bestRespTh);
    
        TGCanvas c = new TGCanvas();
        c.setTitle("Metrics vs Response");
        c.draw(gEff).draw(gPur, "same");
        c.region().showLegend(0.25, 0.25);
        //c.region().axisLimitsY(0.8, 1.01);
    
        return bestRespTh;
    }

    
    
    public static void main(String[] args){        

        //to run
        // /open Level3Particle.java
        // /open Level3Tester_SimulationSIDIS.java
        // Level3Tester_SimulationSIDIS.main(new String[]{});

        String dir = "/Users/tyson/data_repo/trigger_data/sims/claspyth_train/";

        Level3Tester_SimulationSIDIS tester = new Level3Tester_SimulationSIDIS();
        tester.load_trackfinder("clas12rgd.network");
        tester.load_cf("cf_el.network");
        tester.load_LtoStripConv("LtoStrip_convTable.csv");
        tester.load_pider("pid_elPim_fromcfpred.network");//_fromcfpred
        
        DataList preds = tester.getPred(dir,98,90,10000);

        tester.PlotVar(preds, 4, 11, "e-", "Track Response","",0,1,100);

        tester.PlotVar(preds, 1, 11, "e-", "P","[GeV]",0,10,100);
        tester.PlotVar(preds, 2, 11, "e-", "Theta","[Deg]",0,50,100);
        tester.PlotVar(preds, 3, 11, "e-", "Phi","[Deg]",-180,180,100);

        tester.PlotResponse(preds, 0, 11, "e-");
        double bestth=tester.findBestThreshold(preds, 0, 11, 0.995);
        tester.plotVarDep(preds, 0, 11, bestth, true, 1, "P","[GeV]",1,9.0,1.0);
        tester.plotVarDep(preds, 0, 11, bestth, true, 2,"Theta","[Deg]",5.0,35.0,5.);
        tester.plotVarDep(preds, 0, 11, bestth, true, 3,"Phi","[Deg]",-180,180,10.);

    }
}
