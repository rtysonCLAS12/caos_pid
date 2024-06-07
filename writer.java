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
import java.io.PrintWriter;

import j4ml.data.*;
import j4ml.deepnetts.*;
import j4ml.ejml.EJMLModel;
import j4ml.ejml.EJMLModel.ModelType;

/**
 *
 * @author tyson
 */
public class writer {

  EJMLModel cf;
  EJMLModel htccer;
  DataList LtoSconv;
  EJMLModel pider;
  Level3PIDUtils utils = new Level3PIDUtils();

  public writer(){

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

  public void writeOut(String file, String out, int nEvs, Boolean isOut){

    int nFound = 0;
    //use -1 to say we want to read all the input file
    if(nEvs==-1){nFound=-5;}

    // String file=dir;

    HipoReader r = new HipoReader(file);
    HipoWriter w = HipoWriter.create(out, r);

    Event e = new Event();
    Event e_out= new Event();

    // r.getSchemaFactory().show();
    Bank[] banks_nAI = r.getBanks("DC::tdc", "ECAL::adc", "RUN::config", "FTOF::adc", "HTCC::adc",
        "REC::Particle", "REC::Calorimeter", "REC::Cherenkov", "REC::Track", "HitBasedTrkg::HBClusters",
        "HitBasedTrkg::HBTracks", "ECAL::clusters","REC::Scintillator","RUN::config","HTCC::rec","REC::Traj");
    
    
    Bank[] banks = r.getBanks("DC::tdc", "ECAL::adc", "RUN::config", "FTOF::adc", "HTCC::adc",
        "RECAI::Particle", "RECAI::Calorimeter", "RECAI::Cherenkov", "RECAI::Track", "HitBasedTrkg::HBClusters",
        "HitBasedTrkg::HBTracks", "ECAL::clusters","RECAI::Scintillator","RUN::config","HTCC::rec","RECAI::Traj");
    

    CompositeNode tr = new CompositeNode(3210, 2, "i", 1200);
    CompositeNode pt = new CompositeNode(3210, 3, "i", 1200);

     // output
     CompositeNode bcf_out = new CompositeNode(32100, 11, "sf9f9f", 48);
     CompositeNode bcf_out_test = new CompositeNode(32100, 21, "sf9f9f", 48);
     CompositeNode bcf_out_test_nAI = new CompositeNode(32100, 31, "sf9f9f", 48);
     CompositeNode pt_out = new CompositeNode(32100, 3, "sssifffffffsffffffffffffi", 1200);
     CompositeNode pt_out_test = new CompositeNode(32100, 22, "sssifffffffsfffffi", 1200);
     CompositeNode pt_out_test_nAI = new CompositeNode(32100, 32, "sssifffffffsfffffi", 1200);

    while (r.hasNext() && nFound < nEvs) {

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
          //&& cf_pred[9]>0 && sumHTCCADC>0 && !hasAllWrongECALPred(cf_strips) && vars_pid[0]>0
          if (pid_pred[0] > 0.1 && cf_pred[9]>0 && sumHTCCADC>0 && !hasAllWrongECALPred(cf_strips) && vars_pid[0]>0) { 
            outpid = 11;
          } else {
            outpid = -211;
          }
        }

        bcf_out.setRows(i + 1);
        bcf_out.putShort(0, i, (short) i);
        // Predicted ECAL Position
        for (int j = 0; j < 9; j++) {
          bcf_out.putFloat(j + 1, i, cf_strips[j]);
        }
        // Predicted FTOF Position
        bcf_out.putFloat(10, i, cf_pred[9] * 62);
        for (int j = 0; j < 9; j++) {
          //System.out.printf("local pos %f\n",cf_pred[j]);
          bcf_out.putFloat(j + 11, i, cf_pred[j]);
        }

        pt_out.setRows(i + 1);
        pt_out.putShort(0, i, (short) i);
        // Order of Sector and Charge inverted in pt bank compared to tr bank
        pt_out.putShort(1, i, (short) Charge);
        pt_out.putShort(2, i, (short) Sector);
        pt_out.putInt(3, i, outpid);
        // rest of bank is the same
        for (int j = 4; j < 10; j++) {
          pt_out.putFloat(j, i, (float) pt.getDouble(j, i));
        }
        pt_out.putFloat(10, i, resp); // put resp when testing to make plots etc
        pt_out.putShort(11, i, (short) pt.getInt(11, i));
        pt_out.putFloat(12, i, vars_pid[0]);
        pt_out.putFloat(13, i, vars_pid[1]);
        pt_out.putFloat(14, i, vars_pid[2]);
        pt_out.putFloat(15, i, sumHTCCADC); //sumHTCCADC  htcc_pred[0]
        
        //when not using HTCC pred
        pt_out.putFloat(16, i, vars_pid[27]);
        pt_out.putFloat(17, i, vars_pid[28]);
        pt_out.putFloat(18, i, vars_pid[29]);
        pt_out.putFloat(19, i, vars_pid[30]);
        pt_out.putFloat(20, i, vars_pid[31]);
        pt_out.putFloat(21, i, vars_pid[32]);
        pt_out.putFloat(22, i, vars_pid[33]);
        pt_out.putFloat(23, i, vars_pid[34]);
        pt_out.putInt(24, i, banks[13].getInt("event", 0));

        //when using HTCC pred
        /*pt_out.putFloat(16, i, forHTCC_pred[14]);
        pt_out.putFloat(17, i, forHTCC_pred[15]);
        pt_out.putFloat(18, i, forHTCC_pred[16]);
        pt_out.putFloat(19, i, forHTCC_pred[15]);
        pt_out.putFloat(20, i, forHTCC_pred[16]);
        pt_out.putFloat(21, i, forHTCC_pred[17]);
        pt_out.putFloat(22, i, forHTCC_pred[18]);
        pt_out.putFloat(23, i, forHTCC_pred[19]);*/
        

        

        if(Sector==1 && outpid==11){
          nInstaEl++;
        }
        n++;
      }

      int nrows = 0, nRECel=0;
      for (int i = 0; i <  banks[5].getRows(); i++) {

        int[] cf_out_test = new int[10];
        float[] apt_out_test = new float[12];
        float[] ls_test = new float[]{ -2, -2, -2, -2, -2, -2, -2, -2, -2 };
        double track_chi2 = utils.readTestBanks(i, banks[5], banks[6], banks[12], banks[8], cf_out_test,
            apt_out_test, LtoSconv,ls_test);
        float[] extra_out_test = new float[4];
        utils.readTestBanksExtras(i, banks[5], banks[6], banks[12], banks[8], banks[7], extra_out_test);

        if (Math.abs(apt_out_test[11]) >= 2000 && Math.abs(apt_out_test[11]) < 4000 ) {

          bcf_out_test.setRows(nrows + 1);
          bcf_out_test.putShort(0,nrows, (short) i);
          // Predicted ECAL Position
          for (int j = 0; j < 9; j++) {
            bcf_out_test.putFloat(j + 1,nrows, cf_out_test[j]);
          }
          // Predicted FTOF Position
          bcf_out_test.putFloat(10,nrows, cf_out_test[9]);
          for (int j = 0; j < 9; j++) {
            //System.out.printf("local pos %f\n",ls_test[j]);
            bcf_out_test.putFloat(j + 11,nrows, ls_test[j]);
          }

          pt_out_test.setRows(nrows + 1);
          for (int j = 0; j < 3; j++) {
            pt_out_test.putShort(j,nrows, (short) apt_out_test[j]);
          }
          pt_out_test.putInt(3,nrows, (int) apt_out_test[3]);
          // rest of bank is the same
          for (int j = 4; j < 11; j++) {
            pt_out_test.putFloat(j,nrows, apt_out_test[j]);
          }
          pt_out_test.putShort(11,nrows, (short) apt_out_test[11]);
          pt_out_test.putFloat(12,nrows, (float) track_chi2);
          pt_out_test.putFloat(13,nrows, extra_out_test[0]);
          pt_out_test.putFloat(14,nrows, extra_out_test[1]);
          pt_out_test.putFloat(15,nrows, extra_out_test[2]);
          pt_out_test.putFloat(16,nrows, extra_out_test[3]);
          pt_out_test.putInt(17,nrows,banks[13].getInt("event", 0));
          nrows++;

          if(apt_out_test[3]==11 && apt_out_test[2]==1){
            nRECel++;
          }

        }
      }

      int nrows_nAI = 0, nRECel_nAI=0;
      for (int i = 0; i <  banks_nAI[5].getRows(); i++) {

        int[] cf_out_test = new int[10];
        float[] apt_out_test = new float[12];
        float[] ls_test = new float[]{ -2, -2, -2, -2, -2, -2, -2, -2, -2 };
        double track_chi2 = utils.readTestBanks(i, banks_nAI[5], banks_nAI[6], banks_nAI[12], banks_nAI[8], cf_out_test,
            apt_out_test, LtoSconv,ls_test);
        float[] extra_out_test = new float[4];
        utils.readTestBanksExtras(i, banks_nAI[5], banks_nAI[6], banks_nAI[12], banks_nAI[8], banks_nAI[7], extra_out_test);


        if (Math.abs(apt_out_test[11]) >= 2000 && Math.abs(apt_out_test[11]) < 4000 ) {

          bcf_out_test_nAI.setRows(nrows_nAI + 1);
          bcf_out_test_nAI.putShort(0,nrows_nAI, (short) i);
          // Predicted ECAL Position
          for (int j = 0; j < 9; j++) {
            bcf_out_test_nAI.putFloat(j + 1,nrows_nAI, cf_out_test[j]);
          }
          // Predicted FTOF Position
          bcf_out_test_nAI.putFloat(10,nrows_nAI, cf_out_test[9]);
          for (int j = 0; j < 9; j++) {
            bcf_out_test_nAI.putFloat(j + 11,nrows_nAI, ls_test[j]);
          }

          pt_out_test_nAI.setRows(nrows_nAI + 1);
          for (int j = 0; j < 3; j++) {
            pt_out_test_nAI.putShort(j,nrows_nAI, (short) apt_out_test[j]);
          }
          pt_out_test_nAI.putInt(3,nrows_nAI, (int) apt_out_test[3]);
          // rest of bank is the same
          for (int j = 4; j < 11; j++) {
            pt_out_test_nAI.putFloat(j,nrows_nAI, apt_out_test[j]);
          }
          pt_out_test_nAI.putShort(11,nrows_nAI, (short) apt_out_test[11]);
          pt_out_test_nAI.putFloat(12,nrows_nAI, (float) track_chi2);
          pt_out_test_nAI.putFloat(13,nrows_nAI, extra_out_test[0]);
          pt_out_test_nAI.putFloat(14,nrows_nAI, extra_out_test[1]);
          pt_out_test_nAI.putFloat(15,nrows_nAI, extra_out_test[2]);
          pt_out_test_nAI.putFloat(16,nrows_nAI, extra_out_test[3]);
          pt_out_test_nAI.putInt(17,nrows_nAI,banks_nAI[13].getInt("event", 0));
          nrows_nAI++;

          if(apt_out_test[3]==11 && apt_out_test[2]==1){
            nRECel_nAI++;
          }

        }
      }

      int nAll=n+nrows+nrows_nAI;
      if(nAll!=0 && nRECel==0 && nInstaEl>0){

        // pt_out.show();
        // pt_out_test.show();

        e_out.reset();
        e_out.write(pt_out);
        e_out.write(tr);

        e_out.write(bcf_out); // 32100,11

        e_out.write(bcf_out_test);
        e_out.write(pt_out_test);

        e_out.write(bcf_out_test_nAI);
        e_out.write(pt_out_test_nAI);

        for (Bank b : banks) {
          e_out.write(b);
        }

        w.addEvent(e_out);

        // use -1 to say we want to read all the input file
        if (nEvs != -1) {
          nFound++;
        }
       
      }

    }

    System.out.printf("for out "+out+" found %d events\n",nFound);

    w.close();
  }


  public static void main(String[] args) {

    // to run, in order
    // /open Level3Particle.java
    // /open Level3Candidate.java
    // /open Level3Tester.java
    // Level3Tester.main(new String[]{});

    String file = "/Users/tyson/data_repo/trigger_data/rgd/018326/run_18326_1_wAIBanks.h5";
    Boolean outbend=false;

    //String file = "/Users/tyson/data_repo/trigger_data/rgd/018777/run_18777_1_wAIBanks.h5";
    //String file = "/Users/tyson/data_repo/trigger_data/rgd/018777/sorted_out_skim_018777_noPID_wAIBanks.h5";
    //Boolean outbend=true;

    String field="";
    if(outbend){field="_outbending";}

    writer Writer = new writer();
    Writer.load_cf("cf_el_wFTOF"+field+".network");
    Writer.load_LtoStripConv("LtoStrip_convTable.csv");
    Writer.load_pider("old_networks/pid_elNegBG_fromcfpred"+field+".network");// _fromcfpred old_networks/ _wSF _wHTCCPred_wSF old_networks/ old_networks/
    //Writer.load_htccer("htccer_allNeg"+field+".network");

    if(!outbend){field="_inbending";}

    Writer.writeOut(file,"output_test"+field+"_large_v2.h5", -1,outbend); //noRecEl_InstaEl _newNetwork _newNetwork
    //Writer.writeOut(file,"output_noRecEl_InstaEl"+field+"_RfRecooked.h5", 100000,outbend); //_newNetwork


    
  }

  

}