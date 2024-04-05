
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
import j4np.neural.classifier.NeuralClassifierModel;

/**
 *
 * @author tyson
 */
public class Level3Tester_postProcessor {
  NeuralClassifierModel trackfinder = new NeuralClassifierModel();
  EJMLModel cf;
  DataList LtoSconv;
  EJMLModel pider;

  List<float[]> negTracks = new ArrayList<>();
  List<float[]> posTracks = new ArrayList<>();
  List<float[]> all_permutations = new ArrayList<>();

  public Level3Tester_postProcessor() {

  }


  public DataList getPred(String file) {

    int rEv = 0;
    DataList pred = new DataList();

    // String file=dir;

    HipoReader r = new HipoReader(file);

    Event e = new Event();

    CompositeNode tr = new CompositeNode(32100,2,"i",1200);
    CompositeNode pt = new CompositeNode(3210, 3, "i", 1200);
    CompositeNode cf = new CompositeNode(3210, 11, "i", 1200);

    CompositeNode pt_test = new CompositeNode(3210, 22, "i", 1200);
    CompositeNode cf_test = new CompositeNode(3210, 21, "i", 1200);

    while (r.hasNext() ) { //&& rEv<50000
      
      r.next(e);
      e.read(tr,32100,2);
      e.read(pt, 32100, 3);
      e.read(cf, 32100, 11);
      e.read(pt_test, 32100, 22);
      e.read(cf_test, 32100, 21);


      if(pt.getRows()>0 && pt_test.getRows()>0){
        rEv++;
        int hasRecEl=0;
        for(int i=0;i<pt_test.getRows();i++){
          //System.out.printf("vz %f pid %d \n",pt_test.getDouble(9,i),pt_test.getInt(3,i));
          if(pt_test.getDouble(9,i)>-13 && pt_test.getDouble(9,i)<12 && pt_test.getInt(3,i)==11){
            hasRecEl=1;
          }
        }

        double maxPred=0;
        for(int i=0;i<pt.getRows();i++){
          if(pt.getInt(1,i)==-1){
            if(pt.getDouble(10,i)>maxPred){
              maxPred=pt.getDouble(10,i);
            }
          }

        }

        DataEntry de = new DataEntry(
          new double[] { hasRecEl}, // L1trigger[0+sect]
          new double[] { maxPred});
        pred.add(de);

      }
    
      
    }
    return pred;

  }

  public void PlotResponse(DataList pred, int elClass, int desiredPID, String part) {
    int NEvents = pred.getList().size();

    H1F hRespPos = new H1F(part + " in Event", 100, 0, 1);
    hRespPos.attr().setLineColor(2);
    hRespPos.attr().setFillColor(2);
    hRespPos.attr().setLineWidth(3);
    hRespPos.attr().setTitleX("Response");
    H1F hRespNeg = new H1F("No " + part + " in Event", 100, 0, 1);
    hRespNeg.attr().setLineColor(5);
    hRespNeg.attr().setLineWidth(3);
    hRespNeg.attr().setTitleX("Response");
    // Sort predictions into those made on the positive/or negative samples
    for (int i = 0; i < NEvents; i += 1) {
      float[] vars = pred.getList().get(i).floatFirst();
      float[] resp = pred.getList().get(i).floatSecond();

      if (vars[0] == desiredPID) {
        hRespPos.fill(resp[elClass]);
      } else {
        hRespNeg.fill(resp[elClass]);
      }
    }

    TGCanvas c = new TGCanvas();

    c.setTitle("Response");
    c.draw(hRespPos).draw(hRespNeg, "same");
    c.region().showLegend(0.05, 0.95);

  }// End of PlotResponse


  public StatNumber[] getMetrics(DataList pred, int elClass, int desiredPID, double thresh) {

    int nEvents = pred.getList().size();

    double TP = 0, FP = 0, FN = 0;
    for (int i = 0; i < nEvents; i++) {
      float[] vars = pred.getList().get(i).floatFirst();
      float[] resp = pred.getList().get(i).floatSecond();
      if (vars[0] == desiredPID) {
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

    /*
     * System.out.printf("Theres %d electrons in sample\n", nEls);
     * System.out.printf("L1 trigger fired %d times in sample\n", nTrig);
     */
    return calcMetrics(TP, FP, FN);
  }

  public StatNumber[] calcMetrics(double TP, double FP, double FN) {
    StatNumber Pur = new StatNumber(TP, Math.sqrt(TP));
    StatNumber Eff = new StatNumber(TP, Math.sqrt(TP));
    StatNumber FPs = new StatNumber(FP, Math.sqrt(FP));
    StatNumber FNs = new StatNumber(FN, Math.sqrt(FN));
    StatNumber denomPur = new StatNumber(TP, Math.sqrt(TP));
    StatNumber denomEff = new StatNumber(TP, Math.sqrt(TP));
    denomPur.add(FPs);
    denomEff.add(FNs);
    Pur.divide(denomPur);
    Eff.divide(denomEff);
    StatNumber[] mets = new StatNumber[2];
    mets[0] = Pur;
    mets[1] = Eff;
    return mets;
  }



  public double findBestThreshold(DataList pred, int elClass, int desiredPID, double effLow) {

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
    double bestPuratEffLow = 0;
    double bestPurErratEffLow = 0;
    double bestEffatEffLow = 0;
    double bestEffErratEffLow = 0;

    // Loop over threshold on the response
    for (double RespTh = 0.01; RespTh < 0.99; RespTh += 0.01) {
      StatNumber[] metrics = getMetrics(pred, elClass, desiredPID, RespTh);
      double Pur = metrics[0].number();
      double Eff = metrics[1].number();
      double PurErr = metrics[0].error();
      double EffErr = metrics[1].error();
      gPur.addPoint(RespTh, Pur, 0, PurErr);
      gEff.addPoint(RespTh, Eff, 0, EffErr);
      if (Eff > effLow) {
        if (Pur > bestPuratEffLow) {
          bestPuratEffLow = Pur;
          bestPurErratEffLow = PurErr;
          bestEffatEffLow = Eff;
          bestEffErratEffLow = EffErr;
          bestRespTh = RespTh;
        }
      }
    } // Increment threshold on response

    System.out.format(
        "%n Purity %.3f +/- %.4f at Efficiency %.3f +/- %.4f at a threshold on the response of %.3f %n%n",
        bestPuratEffLow, bestPurErratEffLow,bestEffatEffLow,bestEffErratEffLow, bestRespTh);

    TGCanvas c = new TGCanvas();
    c.setTitle("Metrics vs Response");
    c.draw(gEff).draw(gPur, "same");
    c.region().showLegend(0.25, 0.25);
    // c.region().axisLimitsY(0.8, 1.01);

    return bestRespTh;
  }

  public static void main(String[] args) {

    // to run, in order
    // /open Level3Tester_postProcessor.java
    // Level3Tester_postProcessor.main(new String[]{});

    String file = "output_test.h5";
    //String file = "/Users/tyson/data_repo/trigger_data/rgd/018777/run_018777_wAI.h5";
    //String file = "/Users/tyson/data_repo/trigger_data/rgd/018442/run_018442_1_wAI.h5";

    Level3Tester_postProcessor tester = new Level3Tester_postProcessor();

    DataList preds = tester.getPred(file);
    tester.PlotResponse(preds, 0, 1, "e-");
    double bestth = tester.findBestThreshold(preds, 0, 1, 0.99); //0.995

  }

}