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
public class plotter_tracks {

  public plotter_tracks() {

  }

  public int[] hasMatchedAndInd(CompositeNode pt_test, double[] part, double[] matching_part, int charge){

    double Pz_o=part[2];
    int[] matched=new int[2];
    matched[0]=-1;
    matched[1]=-1;

    double bestRes=999;
    int bestInd=-1;

    for(int i = 0; i < pt_test.getRows(); i++){
      if(pt_test.getInt(2,i)==1){
        double resPz=Math.abs(pt_test.getDouble(6,i)-Pz_o)/Math.abs(Pz_o);
        if(pt_test.getInt(1,i)==charge && resPz<0.2){
          if(resPz<bestRes){
            bestRes=resPz;
            bestInd=i;
          }
        }
      }
    }

    if(bestInd!=-1){
      matched[0]=pt_test.getInt(3,bestInd);
      matched[1]=bestInd;
      matching_part[0]=pt_test.getDouble(4,bestInd);
      matching_part[1]=pt_test.getDouble(5,bestInd);
      matching_part[2]=pt_test.getDouble(6,bestInd);
    } else{
      matching_part[0]=0;
      matching_part[1]=0;
      matching_part[2]=0;
    }

    return matched;
  }

  public Boolean passFid(CompositeNode cf, CompositeNode pt,int row,Boolean isREC){
    Boolean pass = true;
    for(int k = 11; k < 14; k++){
      int det =13;
      if(k<4){det=13;}
      else if (k>3 && k <7){det=14;}
      else{det=15;}
      if(isREC){det=det+1;}
      //System.out.printf("k %d & det %d\n",k,det);
      if(cf.getDouble(k,row)<9 && cf.getDouble(k,row)>0 ){ //&& pt.getDouble(det, row)>0
        pass=false;
      }
    }
    return pass;

  }


  public GraphErrors[] plot(String file, int nEv){

    //don't require matching
    //already require that have tracks in given sector difference
    //in rec, this is essentially like matching, but event by event
    //instead of track by track
    Boolean requireMatching=false;

    String endTitle="(Matched Tracks)";
    if(!requireMatching){endTitle="(Matched & Not Matched Tracks)";}

    GraphErrors gEff = new GraphErrors();
    gEff.attr().setMarkerColor(4);
    gEff.attr().setMarkerSize(16);
    gEff.attr().setLineWidth(2);
    gEff.attr().setTitle("Efficiency "+endTitle);
    gEff.attr().setTitleX("Sector difference");
    gEff.attr().setTitleY("Efficiency");

    HipoReader r = new HipoReader(file);

    CompositeNode pt = new CompositeNode(32100,3,"sssifffffffsffffffffffff",1200);
    CompositeNode pt_test = new CompositeNode(32100,22,"sssifffffffsfffff",1200);
    CompositeNode cf = new CompositeNode(32100,11,"i",48);
    CompositeNode cf_test = new CompositeNode(32100,21,"i",48);
      

    Event event = new Event();
    int counter = 0;

    double[] InstaSecDif= new double[4];
    double[] RecSecDif= new double[4];
    for (int l=0; l<4;l++){
      InstaSecDif[l]=0;
      RecSecDif[l]=0;
    }

    while(r.hasNext() && counter<nEv){
      counter++;
      r.next(event);

      //System.out.println("reading");

      event.read(pt,32100,3);
      event.read(pt_test,32100,22);
      event.read(cf,32100,11);
      event.read(cf_test,32100,21);

      int nrows = pt.getRows();
      int nrows_test=pt_test.getRows();
      //int nrowstest_nAI=pt_test_nAI.getRows();

      //System.out.printf(" nrows insta %d nrows ai rec %d \n",nrows,nrowstest);

      double[] InstaSecDif_pEv= new double[4];
      double[] RecSecDif_pEv= new double[4];
      for (int l=0; l<4;l++){
        InstaSecDif_pEv[l]=0;
        RecSecDif_pEv[l]=0;
      }

      for (int i = 0; i < nrows; i++) {
        double Px1=pt.getDouble(4,i);
        double Py1=pt.getDouble(5,i);
        double Pz1=pt.getDouble(6,i);

        int isMatched1=1;
        if(requireMatching){
          double[] el= new double[3];
          el[0]=Px1;
          el[1]=Py1;
          el[2]=Pz1;
          double[] match= new double[3];
          int[] isMatchedArr=hasMatchedAndInd(pt_test, el, match, pt.getInt(1,i));
          isMatched1=isMatchedArr[0];
        }

        if(isMatched1!=-1 && pt.getInt(1,i)!=0 && pt.getInt(3,i)!=11){

          for (int k = i + 1; k < nrows; k++) {

            double Px2=pt.getDouble(4,k);
            double Py2=pt.getDouble(5,k);
            double Pz2=pt.getDouble(6,k);
  
            int isMatched2=1;
            if(requireMatching){
              double[] el= new double[3];
              el[0]=Px2;
              el[1]=Py2;
              el[2]=Pz2;
              double[] match= new double[3];
              int[] isMatchedArr=hasMatchedAndInd(pt_test, el, match, pt.getInt(1,k));
              isMatched2=isMatchedArr[0];
            }
            if(isMatched2!=-1 && pt.getInt(1,k)!=0 && pt.getInt(3,k)!=11 ){
              int sDif=(Math.abs(pt.getInt(2,i)-pt.getInt(2,k)));
              if(sDif==4){sDif=2;};
              if(sDif==5){sDif=1;};
              InstaSecDif_pEv[sDif]++;
            }//check 2nd track is matched & isn't neutral
          }//loop over second part
        }//check 1st track is matched & isn't neutral
      }//loop over 1st track

      for (int i = 0; i < nrows_test; i++) {
        double Px1=pt_test.getDouble(4,i);
        double Py1=pt_test.getDouble(5,i);
        double Pz1=pt_test.getDouble(6,i);

        int isMatched1=1;
        if(requireMatching){
          double[] el= new double[3];
          el[0]=Px1;
          el[1]=Py1;
          el[2]=Pz1;
          double[] match= new double[3];
          int[] isMatchedArr=hasMatchedAndInd(pt, el, match, pt_test.getInt(1,i));
          isMatched1=isMatchedArr[0];
        }

        //
        if(isMatched1!=-1 && pt_test.getInt(1,i)!=0 && pt_test.getInt(3,i)!=11 && pt_test.getDouble(12,i)<350 && Math.abs(pt_test.getDouble(9,i))<20){

          for (int k = i + 1; k < nrows_test; k++) {

            double Px2=pt_test.getDouble(4,k);
            double Py2=pt_test.getDouble(5,k);
            double Pz2=pt_test.getDouble(6,k);
  
            int isMatched2=1;
            if(requireMatching){
              double[] el= new double[3];
              el[0]=Px2;
              el[1]=Py2;
              el[2]=Pz2;
              double[] match= new double[3];
              int[] isMatchedArr=hasMatchedAndInd(pt, el, match, pt_test.getInt(1,k));
              isMatched2=isMatchedArr[0];
            }

            //
            if(isMatched2!=-1 &&  pt_test.getInt(1,k)!=0 && pt_test.getInt(3,k)!=11 && pt_test.getDouble(12,k)<350 && Math.abs(pt_test.getDouble(9,k))<20){
              int sDif=(Math.abs(pt_test.getInt(2,i)-pt_test.getInt(2,k)));
              if(sDif==4){sDif=2;};
              if(sDif==5){sDif=1;};
              RecSecDif_pEv[sDif]++;
            }//check 2nd track is matched & isn't neutral
          }//loop over second part
        }//check 1st track is matched & isn't neutral
      }//loop over 1st track

      for (int l=0; l<4;l++){
        if(RecSecDif_pEv[l]>0){
          RecSecDif[l]+=RecSecDif_pEv[l];
          if(InstaSecDif_pEv[l]>0){
            InstaSecDif[l]+=InstaSecDif_pEv[l];
          }
        }
      }


    }//events

    

    for (int l=0; l<4;l++){

      //can get eff >1 due to more than one Insta electron per event with REC electron
      //don't really want to plot this here because we're not wanting to measure
      //increase in electron efficiency, we want to measure how many times we find
      // an electron when REC has an electron
      if(InstaSecDif[l]>RecSecDif[l]){InstaSecDif[l]=RecSecDif[l];}


      //StatNumber Pur = new StatNumber(hasInstaEl_hasAIEl, Math.sqrt(hasInstaEl_hasAIEl));
      StatNumber Eff = new StatNumber(InstaSecDif[l], Math.sqrt(InstaSecDif[l]));
      //StatNumber denomPur = new StatNumber(hasInstaEl, Math.sqrt(hasInstaEl));
      StatNumber denomEff = new StatNumber(RecSecDif[l], Math.sqrt(RecSecDif[l]));
      //Pur.divide(denomPur);
      Eff.divide(denomEff);
    
      //gPur.addPoint(l Pur.number(), 0, Pur.error());
      gEff.addPoint(l, Eff.number(), 0,Eff.error());
    }

    TGCanvas c = new TGCanvas();
    c.setTitle("Efficiency vs Sector Difference");
    c.draw(gEff);//.draw(gPur, "same");
    //c.region().showLegend(0.65,0.25);

    GraphErrors[] gs = new GraphErrors[2];
    gs[0]=gEff;
    //gs[1]=gPur;
    return gs;

  }


  public static void main(String[] args) {

    String file="output_test_inbending_large.h5"; //_newNetwork

    plotter_tracks pter = new plotter_tracks();
    GraphErrors[] gs_p_matched=pter.plot(file,1000000);


  }

}

