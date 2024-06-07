import j4np.data.base.DataFrame;
import j4np.hipo5.data.Bank;
import j4np.hipo5.data.CompositeNode;
import j4np.hipo5.data.Event;
import j4np.hipo5.data.Node;
import j4np.hipo5.io.HipoReader;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import java.util.List;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
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

import j4ml.data.*;
import j4ml.deepnetts.*;
import j4ml.ejml.EJMLModel;
import j4ml.ejml.EJMLModel.ModelType;

/**
 *
 * @author tyson
 */
public class counter {

  EJMLModel cf;
  EJMLModel htccer;
  DataList LtoSconv;
  EJMLModel pider;
  Level3PIDUtils utils = new Level3PIDUtils();

  public counter(){

  }

  public void load_htccer(String path) {
    htccer = new EJMLModel(path, ModelType.TANH_LINEAR);
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

  public Boolean IntarrayContains(int[] arr, int c) {
    for (int value : arr) {
      if (value == c) {
        return true;
      }
    }
    return false;
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

  public Boolean hasAllWrongECALPred(int[] arr) {
    
    if (arr[0] == -2 && arr[1]==-2 && arr[2]==-2) {
      return true;
    }

    if (arr[3] == -2 && arr[4]==-2 && arr[5]==-2) {
      return true;
    }

    if (arr[6] == -2 && arr[7]==-2 && arr[8]==-2) {
      return true;
    }
    
    return false;
  }

  public void count(String file, int nEvs, Boolean isOut,int trigBitNb){

    int nFound = 0;
    //use -1 to say we want to read all the input file
    if(nEvs==-1){nFound=-5;}

    // String file=dir;

    HipoReader r = new HipoReader(file);

    Event e = new Event();

    // r.getSchemaFactory().show();
    Bank[] banks_nAI = r.getBanks("DC::tdc", "ECAL::adc", "RUN::config", "FTOF::adc", "HTCC::adc",
        "REC::Particle", "REC::Calorimeter", "REC::Cherenkov", "REC::Track", "HitBasedTrkg::HBClusters",
        "HitBasedTrkg::HBTracks", "ECAL::clusters","REC::Scintillator","RUN::config","HTCC::rec");
    
    
    Bank[] banks = r.getBanks("DC::tdc", "ECAL::adc", "RUN::config", "FTOF::adc", "HTCC::adc",
        "RECAI::Particle", "RECAI::Calorimeter", "RECAI::Cherenkov", "RECAI::Track", "HitBasedTrkg::HBClusters",
        "HitBasedTrkg::HBTracks", "ECAL::clusters","RECAI::Scintillator","RUN::config","HTCC::rec");
    

    CompositeNode tr = new CompositeNode(3210, 2, "i", 1200);
    CompositeNode pt = new CompositeNode(3210, 3, "i", 1200);

    double nL1evs=0,nInstaevs=0;
    while (r.hasNext() && nFound < nEvs) {
      nFound++;

      r.nextEvent(e);
      e.read(banks);
      e.read(banks_nAI);

      e.read(tr, 32100, 2);
      e.read(pt, 32100, 3);

      int nInstaEl=0, n=0;
      for (int i = 0; i < tr.getRows(); i++) {

        int Sector = tr.getInt(1, i);
        int Charge = tr.getInt(2, i); 
        float[] track_clusters = new float[6];
        if(isOut){Charge=-1*Charge;} // !!!!!!!!!!!!!!!!!!! for outbending data flip charge

        for (int k = 0; k < 6; k++) {
          track_clusters[k] = (float) tr.getDouble(10 + k, i) / 112;
        }

        float[] cf_pred = new float[11];
        cf.getOutput(track_clusters, cf_pred);
        int[] cf_strips = new int[] { -2, -2, -2, -2, -2, -2, -2, -2, -2 };
        utils.convLtoStrip(cf_pred, cf_strips, LtoSconv);

        //old version with HTCC pred and SF
        float[] vars_pid = new float[35];
        float sumHTCCADC = utils.get_varsPID(cf_pred, cf_strips, track_clusters, vars_pid, banks[1], banks[4], Sector);

        //with HTCC pred and SF
        /*float[] forHTCC_pred= new float[30];
        utils.readHTCCBank_forNphePred(track_clusters, forHTCC_pred, banks[4], Sector);
        float[] htcc_pred = new float[1];
        htccer.getOutput(forHTCC_pred, htcc_pred);

        double Px=pt.getDouble(4,i);
        double Py=pt.getDouble(5,i);
        double Pz=pt.getDouble(6,i);
        float P=(float)Math.sqrt(Px*Px+Py*Py+Pz*Pz);

        float[] vars_pid = new float[31];
        utils.get_varsPID_wHTCCPred_wSF(cf_pred, cf_strips, track_clusters, vars_pid, banks[1], htcc_pred[0],P, Sector);*/
        

        int outpid = 211;
        float resp = 0;
        float[] pid_pred = new float[2];
        if (Charge == -1) {
          pider.getOutput(vars_pid, pid_pred);
          resp = pid_pred[0];
          if (pid_pred[0] > 0.1 && cf_pred[9]>0 && sumHTCCADC>0 && !hasAllWrongECALPred(cf_strips) && vars_pid[0]>0) { //&& !IntarrayContains(cf_strips, -2) 
            outpid = 11;
            nInstaEl++;
          } else {
            outpid = -211;
          }
        }
        n++;

      }

      if(nInstaEl>0){
        nInstaevs++;
      }

      long bits = banks[2].getLong("trigger", 0);
      int[] L1trigger = convertL1Trigger(bits);

      if(L1trigger[trigBitNb]==1){
        nL1evs++;
      }

    }

    double ratio=nInstaevs/nL1evs;
    System.out.printf("\n Got %f L1 trigger and %f L3 trigger events\n",nL1evs,nInstaevs);
    System.out.printf("eg ratio %f\n\n", ratio);
  }


  public static void main(String[] args) {

    // to run, in order
    // /open Level3Particle.java
    // /open Level3Candidate.java
    // /open Level3Tester.java
    // Level3Tester.main(new String[]{});

    String file = "/Users/tyson/data_repo/trigger_data/rgd/018326/run_18326_1_wAIBanks.h5";
    Boolean outbend=false;
    int trigBit=0;

    //String file = "/Users/tyson/data_repo/trigger_data/rgd/018777/run_18777_1_wAIBanks.h5";
    //Boolean outbend=true;
    //int trigBit=7;

    String field="";
    if(outbend){field="_outbending";}

    counter counter = new counter();
    counter.load_cf("cf_el_wFTOF"+field+".network");
    counter.load_LtoStripConv("LtoStrip_convTable.csv");
    counter.load_pider("pid_elNegBG_fromcfpred"+field+".network");// _fromcfpred old_networks/ _wSF _wHTCCPred_wSF old_networks/ old_networks/
    //counter.load_htccer("htccer_allNeg"+field+".network");

    counter.count(file, 1000000,outbend,trigBit); //_newNetwork


    
  }

  

}