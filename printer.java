
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

/**
 *
 * @author tyson
 */
public class printer {

  public printer() {

  }

  public int hasMatched(CompositeNode pt_test, double[] part, double[] matching_part, int charge, int sector){

    double Pz_o=part[2];
    int matched=-1;

    double bestRes=999;
    int bestInd=-1;

    for(int i = 0; i < pt_test.getRows(); i++){
      if(pt_test.getInt(2,i)==1){
        double resPz=Math.abs(pt_test.getDouble(6,i)-Pz_o)/Math.abs(Pz_o);
        if(pt_test.getInt(1,i)==charge && resPz<0.2 && pt_test.getInt(2,i)==sector){
          if(resPz<bestRes){
            bestRes=resPz;
            bestInd=i;
          }
        }
      }
    }

    if(bestInd!=-1){
      matched=pt_test.getInt(3,bestInd);
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

  public void print(String file, String out, Boolean requireMatching, Boolean printToScreen,int nEv){
    // Specify the path where you want to save the CSV file
    Path filePath = Paths.get(out);

    // Delete the file if it already exists
    try {
      Files.deleteIfExists(filePath);
    } catch (IOException e) {
      e.printStackTrace();
    }
    
    

    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(filePath))) {

      HipoReader r = new HipoReader(file); //_noRecEl_InstaEl_

      System.out.println("read "+file);


      CompositeNode tr = new CompositeNode(32100,2,"i",1200);
      CompositeNode pt = new CompositeNode(32100,3,"i",1200);
      CompositeNode cf = new CompositeNode(32100,11,"i",1200);
      CompositeNode pt_test = new CompositeNode(32100,22,"i",1200);
      CompositeNode cf_test = new CompositeNode(32100,21,"i",1200);

      Bank htccrec = r.getBank("HTCC::rec");
      Bank htccadc = r.getBank("HTCC::adc");
      Bank rectraj = r.getBank("RECAI::Traj");

      Event event = new Event();
      int counter = 0;
      while(r.hasNext() && counter<nEv){
        r.next(event);
        event.read(tr,32100,2);
        event.read(pt,32100,3);
        event.read(cf,32100,11);
        event.read(pt_test,32100,22);
        event.read(cf_test,32100,21);

        event.read(htccrec);
        event.read(htccadc);
        event.read(rectraj);

        int nrows = tr.getRows();
        System.out.println("\n--- next event");
        
        int shown=0,isMatched=-1,wrote=0;
        for(int i = 0; i < nrows; i++){
          if(tr.getInt(1,i)==1){
            if(shown==0 && printToScreen){
              System.out.printf("--- InstaRec event %d\n\n",pt.getInt(24,0));
            }
            double Px=pt.getDouble(4,i);
            double Py=pt.getDouble(5,i);
            double Pz=pt.getDouble(6,i);
            double P=Math.sqrt(Px*Px+Py*Py+Pz*Pz);
            
            if(pt.getInt(3,i)==11){
              double[] el= new double[3];
              el[0]=Px;
              el[1]=Py;
              el[2]=Pz;
              double[] match= new double[3];
              isMatched=hasMatched(pt_test, el, match, -1,tr.getInt(1,i));
            }

            //just want unmatched tracks so swap
            if(!requireMatching && pt.getInt(3,i)==11){
              if(isMatched!=-1){
                isMatched=-1;
              } else{
                isMatched=1;  
              }
            }

            if(isMatched!=-1){
              if(printToScreen){
                double Theta = (float) Math.acos((float)Pz / P)*(180.0 / Math.PI);// Math.atan2(Math.sqrt(Px*Px+Py*Py),Pz);
                double Phi = (float) Math.atan2(Py, Px)*(180.0 / Math.PI);
                System.out.printf("sector = %2d, charge = %3d, pid = %3d, resp= %3.3f segments [ ", tr.getInt(1,i),pt.getInt(1,i),pt.getInt(3,i),pt.getDouble(10,i));
                for(int k = 0; k < 6; k++) System.out.printf(" %9.5f ",tr.getDouble(10+k,i));
                System.out.printf("] p %9.5f  %9.5f %9.5f, theta %9.5f, phi %9.5f\n",pt.getDouble(4,i),pt.getDouble(5,i),pt.getDouble(6,i),Theta,Phi);
                System.out.printf("ECAL Strips [ ");
                for(int k = 1; k < 10; k++) System.out.printf("%9.1f ",cf.getDouble(k,i));
                System.out.printf("] FTOF Position %9.1f\npcal ADC = %9.5f, ecin ADC = %9.5f, ecout ADC = %9.5f, sum HTCC ADC %9.5f \n", cf.getDouble(10,i),pt.getDouble(12,i),pt.getDouble(13,i),pt.getDouble(14,i),pt.getDouble(15,i));
                System.out.printf("HTCC PMTS Left: ");
                for(int k = 16; k < 20; k++) System.out.printf(" %9.5f ",pt.getDouble(k,i));
                System.out.printf("\nHTCC PMTS Right: ");
                for(int k = 20; k < 24; k++) System.out.printf(" %9.5f ",pt.getDouble(k,i));
                System.out.println("\n");
                shown++;
              }
              if(wrote==0){
                writer.println(String.valueOf(pt.getInt(24,0)));
                counter++;
                wrote++;
              }
              
            }
          }
          
        }
        
        int nrowstest=pt_test.getRows();
        int shown_r=0;
        for(int i = 0; i < nrowstest; i++){
          if(pt_test.getInt(2,i)==1){
            if(shown_r==0 && printToScreen){
              System.out.printf("--- REC AI event %d\n\n",pt_test.getInt(17,0));
              System.out.println("HTCC adc");
              htccadc.show();
              System.out.println("HTCC rec");
              htccrec.show();
              System.out.println("Traj rec");
              rectraj.show();
            }
            if(isMatched!=-1){
              if(printToScreen){
                double Px=pt_test.getDouble(4,i);
                double Py=pt_test.getDouble(5,i);
                double Pz=pt_test.getDouble(6,i);
                double P=Math.sqrt(Px*Px+Py*Py+Pz*Pz);
                double Theta = (float) Math.acos((float)Pz / P)*(180.0 / Math.PI);// Math.atan2(Math.sqrt(Px*Px+Py*Py),Pz);
                double Phi = (float) Math.atan2(Py, Px)*(180.0 / Math.PI);
                System.out.printf("sector = %2d, charge = %3d, pid = %3d, status = %3d, track chi^2 = %9.5f ", pt_test.getInt(2,i),pt_test.getInt(1,i),pt_test.getInt(3,i),pt_test.getInt(11,i),pt_test.getDouble(12,i));
                System.out.printf("p %9.5f  %9.5f %9.5f, theta %9.5f, phi %9.5f\n",pt_test.getDouble(4,i),pt_test.getDouble(5,i),pt_test.getDouble(6,i),Theta,Phi);
                System.out.printf("ECAL Strips [ ");
                for(int k = 1; k < 10; k++) System.out.printf("%9.1f ",cf_test.getDouble(k,i));
                System.out.printf("] FTOF Position %9.1f \n", cf_test.getDouble(10,i));
                System.out.printf("pcal e = %9.5f, ecin e = %9.5f, ecout e = %9.5f, nb nphe = %9.5f \n\n", pt_test.getDouble(13,i), pt_test.getDouble(14,i), pt_test.getDouble(15,i), pt_test.getDouble(16,i));
                shown_r++;
              }
            }
            
          }
        }
      }
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  public static void main(String[] args) {

    //String field="inbending";
    String field="outbending";
    Boolean useMatched=true;
    Boolean toscreen=true;
    int nEvs=10000;
    if(toscreen){nEvs=100;}


    String uM="_matchedToRECTrack";
    if(useMatched==false){uM="_notMatchedToRECTrack";}

    String file="output_noRecEl_InstaEl_"+field+"_RfRecooked.h5";
    //String out="output_noRec_InstaEl_"+field+"_runList"+uM+".txt";
    String out="bla.txt";

    printer pter = new printer();
    pter.print(file,out,useMatched,toscreen,nEvs);

  }

}