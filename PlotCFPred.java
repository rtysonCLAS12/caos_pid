
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
public class PlotCFPred {
    NeuralClassifierModel trackfinder = new NeuralClassifierModel();
    EJMLModel cf;
    DataList LtoSconv;

    List<float[]> negTracks = new ArrayList<>();
    List<float[]> posTracks = new ArrayList<>();
    List<float[]> all_permutations = new ArrayList<>();

    public PlotCFPred() {

    }

    public void fillHitsHisto(Bank ECALHits, H2F h, int Sector,int det){


        h.attr().setTitleX("Strips");
        h.attr().setTitleY("Layers");

        int l_min=0;
        int l_max=4;
        if(det==1){
            l_min=4;
            l_max=7;            
        } else if(det==2){
            l_min=7;
            l_max=10;            
        }

        for(int k = 0; k < ECALHits.getRows(); k++){
            
            int   sect = ECALHits.getInt("sector", k);
            int  layer = ECALHits.getInt("layer", k);
            int  strip = ECALHits.getInt("component", k);
            double    ADC = ((double) ECALHits.getInt("ADC", k))/5000.;//

            if(ADC>0.0){
                //----------------
                if(sect==Sector){

                    if(det==2){
                        //System.out.printf("sect %d lay %d st %d adc %f \n",sect,layer,strip,ADC);
                    }

                    if(layer>=l_min && layer<l_max){



                        int y  = (layer-1);
                        int x   = strip-1;
                        if(det!=0){
                            y=layer-l_min;
                        }
                        h.fill(x, y,ADC);
                    }
                }
            }
        }
    }

    public void fillStripsHisto(int[] strips, H1F U, H1F V, H1F W,int det,int isTruth){
        int out_start=0;
        if(det==1){
            out_start=3;
        }
        if(det==2){
            out_start=6;
        }
        U.fill(strips[out_start+0]);
        V.fill(strips[out_start+1]);
        W.fill(strips[out_start+2]);

        if(isTruth==1){
            U.attr().setLineColor(2);
            U.attr().setFillColor(2);
            U.attr().setLineWidth(3);
            U.attr().setTitleX("U Cluster Position [strips]");
            V.attr().setLineColor(2);
            V.attr().setFillColor(2);
            V.attr().setLineWidth(3);
            V.attr().setTitleX("V Cluster Position [strips]");
            W.attr().setLineColor(2);
            W.attr().setFillColor(2);
            W.attr().setLineWidth(3);
            W.attr().setTitleX("W Cluster Position [strips]");

        } else{
            U.attr().setLineColor(5);
            U.attr().setLineWidth(3);
            U.attr().setTitleX("U Cluster Position [strips]");
            V.attr().setLineColor(5);
            V.attr().setLineWidth(3);
            V.attr().setTitleX("V Cluster Position [strips]");
            W.attr().setLineColor(5);
            W.attr().setLineWidth(3);
            W.attr().setTitleX("W Cluster Position [strips]");

        }

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

    int[] convLtoStrip(float[] Ls){
        int[] strips = new int[]{-2,-2,-2,-2,-2,-2,-2,-2,-2};
        int[] nStrips = new int[]{68,63,63,37,37,37,37,37,37};
        for(int i=0;i<9;i++){
            Ls[i]=Ls[i]*500;
            float[] conv_det_view = LtoSconv.getList().get(i).floatFirst();
            //System.out.println("conv");
            //System.out.println(Arrays.toString(conv_det_view));
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

    public void plotEx(String file,int nEx){

        // String file=dir;

        HipoReader r = new HipoReader(file);

        Event e = new Event();

        // r.getSchemaFactory().show();
        Bank[] banks = r.getBanks("DC::tdc", "ECAL::adc", "RUN::config", "FTOF::adc", "HTCC::adc",
                "REC::Particle", "REC::Calorimeter", "REC::Cherenkov", "REC::Track", "HitBasedTrkg::HBClusters",
                "HitBasedTrkg::HBTracks", "ECAL::clusters");

        CompositeNode tr = new CompositeNode(3210, 2, "i", 1200);
        CompositeNode pt = new CompositeNode(3210, 3, "i", 1200);

        int rEv=0;
        while (r.hasNext() && rEv<nEx) {
            

            r.nextEvent(e);
            e.read(banks);

            e.read(tr, 32100, 2);
            e.read(pt, 32100, 3);

            // find and initialise particles
        for (int i = 0; i < banks[5].getRows(); i++) {
            Level3Particle part = new Level3Particle();
            part.read_Particle_Bank(i, banks[5]);
            if(part.PIndex!=-1 && part.Charge==-1 && part.PID==11){ //ie part in FD
                part.read_Cal_Bank(banks[6]);
                part.read_HTCC_bank(banks[7]);
                part.find_sector_cal(banks[6]);
                part.find_sector_track(banks[8]);
                part.find_track_clusterIDs(banks[10]);
                part.find_track_clusters(banks[9]);
                part.find_HTCC_ADCs(banks[4]);
                part.set_pid_label();

                //part.print();
                String[] det_names=new String[3];
                det_names[0]="PCAL";
                det_names[1]="ECIN";
                det_names[2]="ECOUT";
                if(part.PCAL_energy!=0 && part.ECIN_energy!=0 && part.ECOUT_energy!=0){
                    for (int det=0;det<3;det++){

                        int nBins=69;
                        double end=68.5;
                        if(det!=0){
                            nBins=37;
                            end=36.5;
                        }

                        float[] cf_pred= new float[9];
                        cf.getOutput(part.track_clusters, cf_pred);
                        int[] cf_strips=convLtoStrip(cf_pred);
                        int[] true_strips=convLtoStrip(part.getNormedCFOut());
                        part.read_cal_bank_from_cf_pred(cf_pred,cf_strips,banks[1]);

                        H1F hUTrue = new H1F("True", nBins, -0.5, end);
                        H1F hUPred = new H1F("Pred", nBins,-0.5,end);
                        H1F hVTrue = new H1F("True", nBins, -0.5, end);
                        H1F hVPred = new H1F("Pred", nBins, -0.5, end);
                        H1F hWTrue = new H1F("True", nBins, -0.5, end);
                        H1F hWPred = new H1F("Pred", nBins, -0.5, end);

                        fillStripsHisto(cf_strips, hUPred, hVPred, hWPred,det,0);
                        fillStripsHisto(true_strips, hUTrue, hVTrue, hWTrue,det,1);

                        H2F hits = new H2F("EC", nBins, 0, end, 3, 0, 3);
                        hits.attr().setTitle(det_names[det]);
                        fillHitsHisto(banks[1], hits, part.Sector,det);

                        TGCanvas c = new TGCanvas();
                        c.view().divide(1, 4);

                        c.setTitle(det_names[det]);
                        c.view().region(0).draw(hits);
                        c.view().region(1).draw(hUTrue).draw(hUPred,"same");
                        c.view().region(2).draw(hVTrue).draw(hVPred,"same");
                        c.view().region(3).draw(hWTrue).draw(hWPred,"same");
                        c.region(1).showLegend(0.05, 0.95);

                    }
                    rEv++;
                }
                
            }
        }

        }

    }

    
    
    public static void main(String[] args){        

        //to run, in order
        // /open Level3Particle.java
        // /open Level3Candidate.java
        // /open PlotCFPred.java
        // PlotCFPred.main(new String[]{});
        String filet = "/Users/tyson/data_repo/trigger_data/rgd/018326/recook_caos_pid/run_018326_4_wAIBanks.h5";

        PlotCFPred conv = new PlotCFPred();
        conv.load_trackfinder("clas12rga.network");
        conv.load_cf("cf_el.network");
        conv.load_LtoStripConv("LtoStrip_convTable.csv");
        
        conv.plotEx(filet,5);//50000

    }
}
