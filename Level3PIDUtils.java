import j4np.hipo5.data.Bank;
import j4np.hipo5.data.CompositeNode;
import j4np.hipo5.data.Event;
import j4np.hipo5.data.Node;
import j4np.hipo5.data.Schema;
import j4np.hipo5.io.HipoReader;
import j4np.hipo5.io.HipoWriter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;
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
import twig.graphics.TGCanvas;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.io.IOException;
import java.io.PrintWriter;

import j4np.data.base.*;
import j4ml.data.*;
import j4ml.deepnetts.*;
import j4ml.ejml.EJMLModel;
import j4ml.ejml.EJMLModel.ModelType;
import j4np.neural.classifier.NeuralClassifierModel;

/**
 *
 * @author tyson
 */
public class Level3PIDUtils {

  public Level3PIDUtils() {
  }

  // convert local U/V/W position in ECAL to strip numbers
  // need strip numbers to read ECAL ADC bank
  void convLtoStrip(float[] Ls, int[] strips, DataList LtoSconv) {

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
  }

  // using ECAL::adc
  public void get_varsPID(float[] cf_pred, int[] cf_strips, float[] track_clusters, float[] vars_for_pid, Bank ECAL_Bank, Bank HTCC_Bank, int sector) {
    
    //Need to initialise some values to 0 as we add to these
    vars_for_pid[0] = 0; //PCAL Energy
    vars_for_pid[1] = 0; //ECIN Energy
    vars_for_pid[2] = 0; //ECOUT Energy
    vars_for_pid[3] = Math.abs(cf_pred[0]) / 500;
    vars_for_pid[4] = Math.abs(cf_pred[1]) / 500;
    vars_for_pid[5] = Math.abs(cf_pred[2]) / 500;
    vars_for_pid[6] = Math.abs(cf_pred[3]) / 500;
    vars_for_pid[7] = Math.abs(cf_pred[4]) / 500;
    vars_for_pid[8] = Math.abs(cf_pred[5]) / 500;
    vars_for_pid[9] = Math.abs(cf_pred[6]) / 500;
    vars_for_pid[10] = Math.abs(cf_pred[7]) / 500;
    vars_for_pid[11] = Math.abs(cf_pred[8]) / 500;
    vars_for_pid[12] = 0 ; //PCAL DU
    vars_for_pid[13] = 0 ; //PCAL DV
    vars_for_pid[14] = 0 ; //PCAL DW
    vars_for_pid[15] = 0 ; //ECIN DU
    vars_for_pid[16] = 0 ; //ECIN DV
    vars_for_pid[17] = 0 ; //ECIN DW
    vars_for_pid[18] = 0 ; //ECOUT DU
    vars_for_pid[19] = 0 ; //ECOUT DV
    vars_for_pid[20] = 0 ; //ECOUT DW

    int n=21;
    //track clusters
    for (float value : track_clusters) {
      vars_for_pid[n] = value;
      n++;
    }

    readECALBank(vars_for_pid,cf_strips,ECAL_Bank,sector);
    //last 8 values are HTCC adc
    readHTCCBank(vars_for_pid,HTCC_Bank,sector);

    /*n=0;
    for (float var : vars_for_pid){
      System.out.printf("var %d val %f, ",n,var);
      n++;
    }
    System.out.println("\n");*/
    
  }

  public void readHTCCBank(float[] vars_for_pid,Bank HTCC_Bank,int sector){
    for(int k = 0; k < HTCC_Bank.getRows(); k++){
      int   sect = HTCC_Bank.getInt("sector", k);
      int  layer = HTCC_Bank.getInt("layer", k); //1 or 2
      int  component = HTCC_Bank.getInt("component", k); //1-4
      float    ADC = HTCC_Bank.getInt("ADC", k);

      //double energy = (ADC / 5000.0); // &&energy<1.0

      if (ADC >= 0.0 && sect == sector) {
          int index = ((layer - 1) * 4 + component) - 1;
          if(index>=0 && index<8){
            vars_for_pid[27+index]=ADC/35000;
          }   
      }
    }
  }

  public void readECALBank(float[] vars_for_pid, int[] cf_strips,Bank ECAL_Bank,int sector){

    //read cal bank
    for (int k = 0; k < ECAL_Bank.getRows(); k++) {
      int sect = ECAL_Bank.getInt("sector", k);
      int layer = ECAL_Bank.getInt("layer", k);
      int strip = ECAL_Bank.getInt("component", k);
      float ADC = ECAL_Bank.getInt("ADC", k);
      if (ADC > 0.0) {
        // ----------------
        if (sect == sector) {
          if (layer == 1) {
            if (strip > (cf_strips[0] - 4) && strip < (cf_strips[0] + 4)) {
              vars_for_pid[12] += 1.0 / 16.0;
              vars_for_pid[0] += ADC / 150000.0;
            }
          } else if (layer == 2) {
            if (strip > (cf_strips[1] - 4) && strip < (cf_strips[1] + 4)) {
              vars_for_pid[13] += 1.0 / 16.0;
              vars_for_pid[0] += ADC / 150000.0;
            }
          } else if (layer == 3) {
            if (strip > (cf_strips[2] - 4) && strip < (cf_strips[2] + 4)) {
              vars_for_pid[14] += 1.0 / 16.0;
              vars_for_pid[0] += ADC / 150000.0;
            }
          } else if (layer == 4) {
            if (strip > (cf_strips[3] - 4) && strip < (cf_strips[3] + 4)) {
              vars_for_pid[15] += 1.0 / 16.0;
              vars_for_pid[1] += ADC / 150000.0;
            }
          } else if (layer == 5) {
            if (strip > (cf_strips[4] - 4) && strip < (cf_strips[4] + 4)) {
              vars_for_pid[16] += 1.0 / 16.0;
              vars_for_pid[1] += ADC / 150000.0;
            }
          } else if (layer == 6) {
            if (strip > (cf_strips[5] - 4) && strip < (cf_strips[5] + 4)) {
              vars_for_pid[17] += 1.0 / 16.0;
              vars_for_pid[1] += ADC / 150000.0;
            }
          } else if (layer == 7) {
            if (strip > (cf_strips[6] - 4) && strip < (cf_strips[6] + 4)) {
              vars_for_pid[18] += 1.0 / 16.0;
              vars_for_pid[2] += ADC / 150000;
            }
          } else if (layer == 8) {
            if (strip > (cf_strips[7] - 4) && strip < (cf_strips[7] + 4)) {
              vars_for_pid[19] += 1.0 / 16.0;
              vars_for_pid[2] += ADC / 150000;
            }
          } else if (layer == 9) {
            if (strip > (cf_strips[8] - 4) && strip < (cf_strips[8] + 4)) {
              vars_for_pid[20] += 1.0 / 16.0;
              vars_for_pid[2] += ADC / 150000;
            }
          }
        }
      }
    }

  }

  public Boolean readTestBanks(int pindex,Bank bpart,Bank bcal, Bank bscint, Bank btrack,int[] cf_out, float[] pt_out,DataList LtoSConv){

    for(int i=0;i<10;i++){
      cf_out[i]=0; //intialise as strip finder sometimes doesn't find strip
    }

    int Track_nSL=6;
    double Track_chi2=0;
                
    for(int i=0;i<btrack.getRows();i++){
      int pind=btrack.getInt("pindex",i);
      if(pind==pindex){
        pt_out[2]=btrack.getInt("sector",i);
        Track_chi2=btrack.getFloat("chi2", i);
        if((btrack.getInt("status", i) & 0b101010101010)!=0){
          Track_nSL=5;
        }
      }
    }

    for(int i=0;i<bscint.getRows();i++){
      int pind=bscint.getInt("pindex",i);
      int layer=bscint.getInt("layer",i);
      if(pind==pindex && layer==2){
        cf_out[9]=bscint.getInt("component",i);
      }
    }

    float[] ls = new float[9];
    for(int i=0; i<bcal.getRows();i++){
      int pind=bcal.getInt("pindex",i);
      int layer=bcal.getInt("layer",i);
      if(pind==pindex){
        ls[layer-1]=bcal.getFloat("lu", i)/500; //L to S conv assumes ls between 0 and 1
        ls[(layer-1)+1]=bcal.getFloat("lv", i)/500;
        ls[(layer-1)+2]=bcal.getFloat("lw", i)/500;
        
      }
    }
    convLtoStrip(ls, cf_out, LtoSConv);

    pt_out[0]=pindex;
    pt_out[1]=bpart.getInt("charge", pindex);
    pt_out[3]=bpart.getInt("pid", pindex);
    pt_out[4]=bpart.getFloat("px", pindex);
    pt_out[5]=bpart.getFloat("py", pindex);
    pt_out[6]=bpart.getFloat("pz", pindex);
    pt_out[7]=bpart.getFloat("vx", pindex);
    pt_out[8]=bpart.getFloat("vy", pindex);
    pt_out[9]=bpart.getFloat("vz", pindex);
    pt_out[10]=bpart.getFloat("vt", pindex);
    pt_out[11]=bpart.getInt("status", pindex);

    Boolean goodTrack=false;
    if(Track_nSL==6 && Track_chi2<350){
      goodTrack=true;
    }
    return goodTrack;

  }

  public static void main(String[] args) {

  }

}