
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
import twig.data.StatNumber;
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

/**
 *
 * @author tyson
 */
public class Level3Histogrammer {
  EJMLModel cf;
  DataList LtoSconv;
  EJMLModel pider;

  List<float[]> negTracks = new ArrayList<>();
  List<float[]> posTracks = new ArrayList<>();
  List<float[]> all_permutations = new ArrayList<>();

  public Level3Histogrammer() {

  }

  public void load_cf(String path) {
    cf = new EJMLModel(path, ModelType.TANH_LINEAR);
  }

  public void load_LtoStripConv(String path) {
    LtoSconv = DataList.fromCSV("LtoStrip_convTable.csv",
        DataList.range(0, 69), DataList.range(0, 1));
  }

  public void load_pider(String path) {
    pider = new EJMLModel(path, ModelType.SOFTMAX);
  }

  public Boolean arrayContains(float[] arr, int c) {
    for (double value : arr) {
      if (value == c) {
        return true;
      }
    }
    return false;
  }

  int[] convLtoStrip(float[] Ls) {
    int[] strips = new int[] { -2, -2, -2, -2, -2, -2, -2, -2, -2 };
    int[] nStrips = new int[] { 68, 63, 63, 37, 37, 37, 37, 37, 37 };
    for (int i = 0; i < 9; i++) {
      Ls[i] = Ls[i] * 500;
      float[] conv_det_view = LtoSconv.getList().get(i).floatFirst();
      for (int str = 0; str < (nStrips[i]); str++) {
        // in PCAL V, W take distance from other end
        // System.out.printf("wtf %d %d \n",i,str);
        if (i == 1 || i == 2) {
          if (Ls[i] < conv_det_view[str] && Ls[i] >= conv_det_view[str + 1]) {
            strips[i] = str + 1; // add 1 as the upper limit applies to strip+1 due to inverted order
          }
        } else {
          if (Ls[i] >= conv_det_view[str] && Ls[i] < conv_det_view[str + 1]) {
            strips[i] = str;
          }
        }
      }
    }
    return strips;
  }

  public Boolean arrayNotEmpty(float[] arr) {
    for (double value : arr) {
      if (value != 0) {
        return true;
      }
    }
    return false;
  }

  public List<Level3Candidate> getCandidatesInSector(List<Level3Candidate> cs, int Sector){
    List<Level3Candidate> scs = new ArrayList<Level3Candidate>();
    for(Level3Candidate c : cs){
      if(c.Sector==Sector){
        scs.add(c);
      }
    }
    return scs;

  }

  public List<Level3Particle> getParticlesInSector(List<Level3Particle> cs, int Sector){
    List<Level3Particle> scs = new ArrayList<Level3Particle>();
    for(Level3Particle c : cs){
      if(c.Sector==Sector){
        scs.add(c);
      }
    }
    return scs;

  }

  public double particlesGetQ2(List<Level3Particle> ps,double beamE){

    int bestInd=-1,bestElInd=-1,i=0;
    double bestP=0,bestElP=0;
    for(Level3Particle p : ps){
      if(p.P>bestP){
        bestP=p.P;
        bestInd=i;
      }
      if(p.Charge == -1 && p.PID==11 && arrayNotEmpty(p.HTCC_adcs) && p.Vz < 12 && p.Vz > -13) {
        if(p.P>bestElP){
          bestElP=p.P;
          bestElInd=i;
        }
      }
      i++;
    }

    if(bestElInd!=-1){
      bestInd=bestElInd;
    }

    double Q2=0;
    if(bestInd!=-1){
      Level3Particle p=ps.get(bestInd);
      double E=Math.sqrt(p.P*p.P+p.getM(11)*p.getM(11));
      Q2 = 2 * beamE * E* (1. - p.Pz/p.P);
    }
    
    return Q2;
  }

  public double candidatesGetQ2(List<Level3Candidate> ps,double beamE,double respTh){

    int bestInd=-1,bestElInd=-1,i=0;
    double bestP=0,bestElP=0;
    for(Level3Candidate p : ps){
      if(p.Pred_P>bestP){
        bestP=p.Pred_P;
        bestInd=i;
      }
      if(p.Pred_Charge == -1 && p.pid_resp > respTh && arrayNotEmpty(p.HTCC_adcs) && p.Vz < 12 && p.Vz > -13) {
        if(p.Pred_P>bestElP){
          bestElP=p.Pred_P;
          bestElInd=i;
        }
      }
      i++;
    }

    if(bestElInd!=-1){
      bestInd=bestElInd;
    }

    double Q2=0;
    if(bestInd!=-1){
      Level3Candidate p=ps.get(bestInd);
      double E=Math.sqrt(p.Pred_P*p.Pred_P+p.getM(11)*p.getM(11));
      Q2 = 2 * beamE * E* (1. - p.Pred_Pz/p.Pred_P);
    }
    
    return Q2;
  }


  // 0 1 2 3 4
  // 5 6 7 8 9
  // 10 11
  // Bank[] banks = r.getBanks("DC::tdc", "ECAL::adc", "RUN::config", "FTOF::adc",
  // "HTCC::adc",
  // "REC::Particle", "REC::Calorimeter",
  // "REC::Cherenkov","REC::Track","HitBasedTrkg::HBClusters",
  // "HitBasedTrkg::HBTracks","ECAL::clusters");
  public List<Level3Candidate> getCandidates(Bank[] banks, CompositeNode tr, CompositeNode pt, Boolean isOut) {

    // Instant start_time = Instant.now();

    List<Level3Candidate> ps = new ArrayList<Level3Candidate>();
    for (int i = 0; i < tr.getRows(); i++) {
      // find and initialise Candidates
      Level3Candidate part = new Level3Candidate();
      // part.setShow(true);
      part.find_AI_tracks(tr, pt, i,isOut);
      part.find_ClosestRECParticle_fromPredP(banks[5], 999, 999, 999);
      part.read_Cal_Bank(banks[6]);
      part.find_HTCC_ADCs(banks[4]);
      part.set_pid_label();

      float[] cf_pred = new float[11];
      cf.getOutput(part.track_clusters, cf_pred);
      int[] cf_strips = convLtoStrip(cf_pred);
      part.read_cal_bank_from_cf_pred(cf_pred, cf_strips, banks[1]);
      float[] pid_pred = new float[2];
      pider.getOutput(part.get_vars_forpid(), pid_pred); //_noTrack
      part.setPidResp(pid_pred[0]);
      // part.print();
      // part.printbanks(tr,pt,banks[5],banks[6],banks[4]);

      ps.add(part);
    }

    // Instant finish_time = Instant.now();
    // timeElapsed += Duration.between(start_time, finish_time).toMillis();

    return ps;
  }

  // 0 1 2 3 4
  // 5 6 7 8 9
  // 10 11
  // Bank[] banks = r.getBanks("DC::tdc", "ECAL::adc", "RUN::config", "FTOF::adc","HTCC::adc",
  // "REC::Particle", "REC::Calorimeter", "REC::Cherenkov","REC::Track","HitBasedTrkg::HBClusters",
  // "HitBasedTrkg::HBTracks","ECAL::clusters");
  public List<Level3Particle> getRECParts(Bank[] banks){

    List<Level3Particle> ps = new ArrayList<Level3Particle>();

    for (int i = 0; i < banks[5].getRows(); i++) {
        Level3Particle part = new Level3Particle();
        part.read_Particle_Bank(i, banks[5]);
        if(part.PIndex!=-1){ //ie part in FD
            part.read_Cal_Bank(banks[6]);
            part.read_HTCC_bank(banks[7]);
            part.find_sector_cal(banks[6]);
            part.find_sector_track(banks[8]);
            part.find_track_clusterIDs(banks[10]);
            part.find_track_clusters(banks[9]);
            part.find_HTCC_ADCs(banks[4]);
            part.set_pid_label();

            float[] cf_pred= new float[11];
            cf.getOutput(part.track_clusters, cf_pred);
            int[] cf_strips=convLtoStrip(cf_pred);
            part.read_cal_bank_from_cf_pred(cf_pred,cf_strips,banks[1]);
            float[] pid_pred =new float[2];
            pider.getOutput(part.get_vars_forpid(), pid_pred);
            part.setPidResp(pid_pred[0]);

            ps.add(part);
        }
    }

    return ps;
  }

  public static int[] convertL1Trigger(long bits) {
    int[] trigger = new int[32];

    // System.out.printf("%X - %X\n", bits,bits&0xF);
    for (int i = 0; i < trigger.length; i++) {
      trigger[i] = 0;
      if (((bits >> i) & (1L)) != 0L)
        trigger[i] = 1;
      // System.out.println(Arrays.toString(trigger));
    }
    return trigger;
  }

  public DataList getPred(String file, int nPart, double respTh, Boolean isOut) {

    int nFound = 0, nL1=0, nNoL1=0;
    DataList pred = new DataList();

    double beamE=10.532;
    double Q2Cut=1.2;

    int L1TrigStart=0;
    if(isOut){L1TrigStart=7;}

    // String file=dir;

    HipoReader r = new HipoReader(file);

    Event e = new Event();

    // r.getSchemaFactory().show();
    Bank[] banks = r.getBanks("DC::tdc", "ECAL::adc", "RUN::config", "FTOF::adc", "HTCC::adc",
        "REC::Particle", "REC::Calorimeter", "REC::Cherenkov", "REC::Track", "HitBasedTrkg::HBClusters",
        "HitBasedTrkg::HBTracks", "ECAL::clusters");
    
    Bank[]  AIbanks = r.getBanks("DC::tdc", "ECAL::adc", "RUN::config", "FTOF::adc", "HTCC::adc",
        "RECAI::Particle", "RECAI::Calorimeter", "RECAI::Cherenkov", "RECAI::Track", "HitBasedTrkg::HBClusters",
        "HitBasedTrkg::HBTracks", "ECAL::clusters");
    

    CompositeNode tr = new CompositeNode(3210, 2, "i", 1200);
    CompositeNode pt = new CompositeNode(3210, 3, "i", 1200);

    while (r.hasNext() && nFound < nPart) {

      r.nextEvent(e);
      e.read(banks);
      e.read(AIbanks);

      e.read(tr, 32100, 2);
      e.read(pt, 32100, 3);

      List<Level3Candidate> Candidates = getCandidates(AIbanks, tr, pt, isOut);
      List<Level3Particle> RECparticles = getRECParts(banks);
      List<Level3Particle> AIRECparticles = getRECParts(AIbanks);

      for (int sect = 1; sect < 7; sect++) {

        List<Level3Candidate> SectorCandidates = getCandidatesInSector(Candidates,sect);
        List<Level3Particle> SectorRECparticles = getParticlesInSector(RECparticles,sect);
        List<Level3Particle> SectorAIRECparticles = getParticlesInSector(AIRECparticles,sect);

        double rTOd= (180.0 / Math.PI);

        long bits = banks[2].getLong("trigger", 0);
        int[] L1trigger = convertL1Trigger(bits);
        if(L1trigger[L1TrigStart + sect]==1){
          nL1++;
        } else {
          nNoL1++;
        }

        double InstaQ2=candidatesGetQ2(SectorCandidates, beamE, respTh);
        double RECQ2 = particlesGetQ2(SectorRECparticles, beamE);
        double RECAIQ2 = particlesGetQ2(SectorAIRECparticles, beamE);

        //want to pass Q2 cut event when no Q2 cut in L1 trigger
        if(!isOut){
          InstaQ2=99999;
          RECAIQ2=99999;
          RECQ2=99999;
        }
        
        for (Level3Candidate p : SectorCandidates) {
          int isEl=0;
          if(p.Pred_Charge == -1 && p.pid_resp > respTh && arrayNotEmpty(p.HTCC_adcs) && p.Vz < 12 && p.Vz > -13) {
            isEl=1;
          }
          double P=p.Pred_P;
          double Theta=p.Pred_Theta*rTOd;
          double Phi=p.Pred_Phi*rTOd;
          double PCALADC=p.PCAL_energy_fcf;
          double TotADC=p.PCAL_energy_fcf+p.ECIN_energy_fcf+p.ECOUT_energy_fcf;
          double HTCCADC=0;
          for (double adc:p.HTCC_adcs){
            HTCCADC+=adc;
          }
          if(InstaQ2>Q2Cut){
            DataEntry de = new DataEntry(
              new double[] {P,Theta,Phi,PCALADC,TotADC,HTCCADC,p.pid_resp,p.PID}, // L1trigger[0+sect]
              new double[] {0,isEl,L1trigger[L1TrigStart + sect]});
            pred.add(de);
          }
          
          
        }

        for (Level3Particle p : SectorRECparticles) {
          int isEl=0;
          if(p.Charge == -1 && p.PID==11 && arrayNotEmpty(p.HTCC_adcs) && p.Vz < 12 && p.Vz > -13) {
            isEl=1;
          }
          double P=p.P;
          double Theta=p.Theta*rTOd;
          double Phi=p.Phi*rTOd;
          double PCALADC=p.PCAL_energy_fcf;
          double TotADC=p.PCAL_energy_fcf+p.ECIN_energy_fcf+p.ECOUT_energy_fcf;
          double HTCCADC=0;
          for (double adc:p.HTCC_adcs){
            HTCCADC+=adc;
          }
          if(RECQ2>Q2Cut){
            DataEntry de = new DataEntry(
              new double[] {P,Theta,Phi,PCALADC,TotADC,HTCCADC,p.pid_resp,p.PID}, // L1trigger[0+sect]
              new double[] {1,isEl,L1trigger[L1TrigStart + sect]});
            pred.add(de);
            nFound++;
          }
          
          
        }

        for (Level3Particle p : SectorAIRECparticles) {
          int isEl=0;
          if(p.Charge == -1 && p.PID==11 && arrayNotEmpty(p.HTCC_adcs) && p.Vz < 12 && p.Vz > -13) {
            isEl=1;
          }
          double P=p.P;
          double Theta=p.Theta*rTOd;
          double Phi=p.Phi*rTOd;
          double PCALADC=p.PCAL_energy_fcf;
          double TotADC=p.PCAL_energy_fcf+p.ECIN_energy_fcf+p.ECOUT_energy_fcf;
          double HTCCADC=0;
          for (double adc:p.HTCC_adcs){
            HTCCADC+=adc;
          }
          if(RECAIQ2>Q2Cut){
            DataEntry de = new DataEntry(
              new double[] {P,Theta,Phi,PCALADC,TotADC,HTCCADC,p.pid_resp,p.PID}, // L1trigger[0+sect]
              new double[] {2,isEl,L1trigger[L1TrigStart + sect]});
            pred.add(de);
          }
          
          
        }

      }

    }

    System.out.printf("\n Nb L1 %d, Nb No L1 %d \n",nL1,nNoL1);

    return pred;

  }

  //whichtoplot: 0 - all, 1 - only InstaRec, 2 - only REC
  //withWithoutTrig: 0 - when events don't pass L1 Trigger, 1 - when events pass L1 Trigger, 2 - both
  public void PlotVar(DataList pred, int varIndex, int whichtoplot, int withWithoutTrig, String varName, String varUnits,
      double low, double high, int nBins) {
    int NEvents = pred.getList().size();

    int nEls=0,nNotEls=0,nRECEls=0,nRECNotEls=0,nRECAIEls=0,nRECAINotEls=0;

    H1F hPosInsta = new H1F("IRec e-", nBins, low, high);
    hPosInsta.attr().setLineColor(2);
    //hPosInsta.attr().setFillColor(2);
    hPosInsta.attr().setLineWidth(3);
    hPosInsta.attr().setTitleX(varName + " " + varUnits);
    H1F hNegInsta = new H1F("IRec !e-", nBins, low, high);
    hNegInsta.attr().setLineColor(5);
    hNegInsta.attr().setLineWidth(3);
    hNegInsta.attr().setTitleX(varName + " " + varUnits);

    H1F hPosREC = new H1F("REC e-", nBins, low, high);
    hPosREC.attr().setLineColor(6);
    //hPosREC.attr().setFillColor(6);
    hPosREC.attr().setLineWidth(3);
    hPosREC.attr().setTitleX(varName + " " + varUnits);
    H1F hNegREC = new H1F("REC !e-", nBins, low, high);
    hNegREC.attr().setLineColor(3);
    hNegREC.attr().setLineWidth(3);
    hNegREC.attr().setTitleX(varName + " " + varUnits);

    H1F hPosRECAI = new H1F("RECAI e-", nBins, low, high);
    hPosRECAI.attr().setLineColor(4);
    //hPosRECAI.attr().setFillColor(4);
    hPosRECAI.attr().setLineWidth(3);
    hPosRECAI.attr().setTitleX(varName + " " + varUnits);
    H1F hNegRECAI = new H1F("RECAI !e-", nBins, low, high);
    hNegRECAI.attr().setLineColor(8);
    hNegRECAI.attr().setLineWidth(3);
    hNegRECAI.attr().setTitleX(varName + " " + varUnits);

    //DataEntry de = new DataEntry(
    //new double[] {P,Theta,Phi,PCALADC,TotADC,HTCCADC,p.pid_resp,p.PID}, // L1trigger[0+sect]
    //new double[] {isInstaOrRECorRECAI,isEl,L1trigger[0 + sect]});

    // Sort predictions into those made on the positive/or negative samples
    for (int i = 0; i < NEvents; i += 1) {
      float[] vars = pred.getList().get(i).floatFirst();
      float[] label = pred.getList().get(i).floatSecond();

      if(withWithoutTrig==2){
        if(label[0]==0){
          if(label[1]==0){
            hNegInsta.fill(vars[varIndex]);
          } else{
            hPosInsta.fill(vars[varIndex]);
          }
        } else if(label[0]==1){
          if(label[1]==0){
            hNegREC.fill(vars[varIndex]);
          } else{
            hPosREC.fill(vars[varIndex]);
          }
        } else if(label[0]==2){
          if(label[1]==0){
            hNegRECAI.fill(vars[varIndex]);
          } else{
            hPosRECAI.fill(vars[varIndex]);
          }
        }
      } else if(withWithoutTrig==label[2]){// 0 or 1
        if(label[0]==0){
          if(label[1]==0){
            hNegInsta.fill(vars[varIndex]);
            nNotEls++;
          } else{
            hPosInsta.fill(vars[varIndex]);
            nEls++;
          }
        } else if(label[0]==1){
          if(label[1]==0){
            hNegREC.fill(vars[varIndex]);
            nRECNotEls++;
          } else{
            hPosREC.fill(vars[varIndex]);
            nRECEls++;
          }
        } else if(label[0]==2){
          if(label[1]==0){
            hNegRECAI.fill(vars[varIndex]);
            nRECAINotEls++;
          } else{
            hPosRECAI.fill(vars[varIndex]);
            nRECAIEls++;
          }
        }

      }
  
    }

    String trigName="With & Without L1 Trigger";
    if(withWithoutTrig==0){
      trigName="Without L1 Trigger";
    } else if(withWithoutTrig==1){
      trigName="With L1 Trigger";
    }

    TGCanvas c = new TGCanvas();
    c.setTitle(varName+" "+trigName);
    if(whichtoplot==0){
      c.draw(hPosInsta).draw(hNegInsta, "same").draw(hPosREC,"same").draw(hNegREC, "same").draw(hPosRECAI,"same").draw(hNegRECAI, "same");
    } else if(whichtoplot==1){
      c.draw(hPosInsta).draw(hNegInsta, "same");
    } if(whichtoplot==2){
      c.draw(hPosREC).draw(hNegREC, "same").draw(hPosRECAI,"same").draw(hNegRECAI, "same");
    }
    
    c.region().showLegend(0.7, 0.95);

    System.out.printf("\n "+trigName+" nb IREC e- %d, nb IREC !e- %d, nb REC e- %d, nb REC !e- %d, nb RECAI e- %d, nb RECAI !e- %d\n",nEls,nNotEls,nRECEls,nRECNotEls,nRECAIEls,nRECAINotEls);

  }// End of PlotResponse

  public static void main(String[] args) {

    // to run, in order
    // /open Level3Particle.java
    // /open Level3Candidate.java
    // /open Level3Histogrammer.java
    // Level3Histogrammer.main(new String[]{});

    String file = "/Users/tyson/data_repo/trigger_data/rgd/018326/run_18326_1_wAIBanks.h5";
    Boolean outbend=false;

    //String file = "/Users/tyson/data_repo/trigger_data/rgd/018777/run_18777_3_wAIBanks.h5";
    //Boolean outbend=true;


    String field="";
    if(outbend){field="_outbending";}

    Level3Histogrammer tester = new Level3Histogrammer();
    tester.load_cf("cf_el_wFTOF"+field+".network");
    tester.load_LtoStripConv("LtoStrip_convTable.csv");
    tester.load_pider("old_networks/pid_elNegBG_fromcfpred"+field+".network");// _fromcfpred
    
    DataList pred = tester.getPred(file,1000000 ,0.06,outbend); 

    //whichtoplot: 0 - all, 1 - only InstaRec, 2 - only REC
    //withWithoutTrig: 0 - when events don't pass L1 Trigger, 1 - when events pass L1 Trigger, 2 - both
    //PlotVar(DataList pred, int varIndex, int whichtoplot, int withWithoutTrig, String varName, String varUnits, double low, double high, int nBins) 
    //pred:
    //new double[] {P,Theta,Phi,PCALADC,TotADC,HTCCADC,p.pid_resp,p.PID}, // L1trigger[0+sect]
    //new double[] {isInstaOrRECorRECAI,isEl,L1trigger[0 + sect]});
    tester.PlotVar(pred, 0, 0, 1, "P", "[GeV]", 0, 10, 100);
    
    tester.PlotVar(pred, 0, 0, 0, "P", "[GeV]", 0, 10, 100);

    //tester.PlotVar(pred, 1, 0, 1, "Theta", "[Deg]", 0, 50, 100);
    tester.PlotVar(pred, 3, 0, 2, "PCAL ADC", "", 0, 150000, 500);
    tester.PlotVar(pred, 4, 0, 2, "Total ECAL ADC", "", 0, 200000, 500);
    tester.PlotVar(pred, 5, 0, 2, "Total HTCC ADC", "", 0, 35000, 500);

    tester.PlotVar(pred, 3, 0, 2, "PCAL ADC (low)", "", 0, 400, 400);
    tester.PlotVar(pred, 4, 0, 2, "Total ECAL ADC (low)", "", 0, 200, 200);
    tester.PlotVar(pred, 5, 0, 2, "Total HTCC ADC (low)", "", 0, 150, 150);

    tester.PlotVar(pred, 5, 0, 2, "Total HTCC ADC (med)", "", 0, 500, 50);

    //tester.PlotVar(pred, 1, 0, 1, "Phi", "[Deg]", -180, 180, 100);

    //tester.PlotVar(pred, 1, 0, 0, "PCAL ADC", "", 0, 200, 200);
    //tester.PlotVar(pred, 1, 0, 0, "Total ECAL ADC", "", 0, 200, 200);
  }
}
