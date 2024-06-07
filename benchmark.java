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
public class benchmark {

  public benchmark() {

  }

  public static void Process(int nThreads){

    String file = "/Users/tyson/data_repo/trigger_data/rgd/018326/run_18326_1_wAIBanks.h5";

    DataStream<HipoReader, HipoWriter, Event> str = new DataStream<HipoReader, HipoWriter, Event>();

    DataFrame<Event> frame = new DataFrame<>();
    HipoReader source = new HipoReader();
    HipoWriter sync = new HipoWriter();

    sync.open("output_1.h5");

    source.open(file);

    DataWorker<HipoReader, Event> worker = new DataWorker<HipoReader, Event>() {

      public final AtomicInteger counter = new AtomicInteger(0);
      Level3PIDUtils utils = new Level3PIDUtils();
      private Schema ECALadc_schema = null;
      private Schema HTCCadc_schema = null;
      EJMLModel cf;
      DataList LtoSconv;
      EJMLModel pider;

      @Override
      public void execute(Event e) {

        //input HTCC ECAL ADC
        Bank bHTCC = new Bank(HTCCadc_schema);
        Bank bECAL = new Bank(ECALadc_schema);
        e.read(bHTCC);
        e.read(bECAL);

        //input track
        CompositeNode tr = new CompositeNode(32100, 2, "i", 1200);
        CompositeNode pt = new CompositeNode(32100, 3, "i", 1200);
        e.read(tr, 32100, 2);
        e.read(pt, 32100, 3);

        //output
        CompositeNode bcf_out = new CompositeNode(32100,11,"sf9f",24);
        CompositeNode pt_out = new CompositeNode(32100,3,"sssifffffffs",1200);

        for (int i = 0; i < tr.getRows(); i++) {

          int Sector=tr.getInt(1,i);
          int Charge=tr.getInt(2,i);
          float[] track_clusters = new float[6];

          for (int k = 0; k < 6; k++) {
            track_clusters[k] = (float) tr.getDouble(10 + k, i) / 112;
          }

          float[] cf_pred = new float[11];
          cf.getOutput(track_clusters, cf_pred);
          int[] cf_strips = new int[] { -2, -2, -2, -2, -2, -2, -2, -2, -2 };
          utils.convLtoStrip(cf_pred,cf_strips,LtoSconv);

          int outpid=211;
          if(Charge==-1){
            float[] vars_pid = new float[35];
            float sumHTCCADC = utils.get_varsPID(cf_pred, cf_strips, track_clusters,vars_pid, bECAL, bHTCC, Sector);
            float[] pid_pred = new float[2];
            pider.getOutput(vars_pid, pid_pred);
            if(pid_pred[0]>0.1 && cf_pred[9]>0 && sumHTCCADC>0  && vars_pid[0]>0 && !hasAllWrongECALPred(cf_strips)){
              outpid=11;
            } else{
              outpid=-211;
            }
          }

          bcf_out.setRows(i+1); 
          bcf_out.putShort( 0, i, (short) i);
          //Predicted ECAL Position
          for(int j=0;j<9;j++){
            bcf_out.putFloat(j+1,i,cf_strips[j]);
          }
          //Predicted FTOF Position
          bcf_out.putFloat(10,i,cf_pred[9]*62);

          pt_out.setRows(i+1); 
          pt_out.putShort(  0, i, (short) i);
          //Order of Sector and Charge inverted in pt bank compared to tr bank
          pt_out.putShort(  1, i, (short) Charge);
          pt_out.putShort(  2, i, (short) Sector);
          pt_out.putInt(  3, i, outpid);
          //rest of bank is the same
          for(int j=4;j<11;j++){
            pt_out.putFloat(  j, i, (float) pt.getDouble(j, i));
          }
          pt_out.putShort( 11, i, (short) pt.getInt(11, i));

        }

        //below fails with: 
        //error replacing node ( 32100,     3) dues to size difference 126,132
        //e.replace(32100, 3, pt_out);

        //Messing with outstream because e.remove prints LENGTH = somelength
        //this was annoying me
        //PrintStream original = System.out;
        //System.setOut(new PrintStream(OutputStream.nullOutputStream()));

        e.remove(32100, 3);

        //System.setOut(original);
        e.write(pt_out);

        e.write(bcf_out); //32100,11
        counter.set(counter.intValue()+1);

        /*try {
          Thread.sleep(25);
        } catch (InterruptedException ex) {
          Logger.getLogger(DataStream.class.getName()).log(Level.SEVERE, null, ex);
        }*/

      }

      @Override
      public boolean init(HipoReader src) {
        // ECAL ADC schema used for prediction
        ECALadc_schema = src.getSchemaFactory().getSchema("ECAL::adc");
        HTCCadc_schema = src.getSchemaFactory().getSchema("HTCC::adc");

        // load networks
        cf = new EJMLModel("cf_el_wFTOF.network", ModelType.TANH_LINEAR);
        LtoSconv = DataList.fromCSV("LtoStrip_convTable.csv",
            DataList.range(0, 69), DataList.range(0, 1));
        pider = new EJMLModel("pid_elNegBG_fromcfpred.network", ModelType.SOFTMAX);
        return true;
      }

      //returns false when all predicted strips are -2
      //eg the track goes out of acceptance of one layer of calorimeter
      //this often indicates a bad track
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

      public void show() {
        System.out.printf(" counter value = %d\n", counter.intValue());
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

    };

    for (int i = 0; i < 8; i++) {
      frame.addEvent(new Event());
    }
    str.threads(nThreads);


    //double timeElapsed=0;
    //Instant start_time = Instant.now();

    str.withSource(source).withOutput(sync).withFrame(frame).consumer(worker).run();
    str.show();

    //Instant finish_time = Instant.now();
    //timeElapsed += Duration.between(start_time, finish_time).getSeconds();//.toMillis();
    //System.out.printf("Took %f s\n",timeElapsed);
  }

  public static int getNEvs(String file){
    HipoReader r = new HipoReader(file);

    CompositeNode pt = new CompositeNode(32100,3,"sssifffffffsffffffffffff",1200);

    Event event = new Event();
    int counter = 0;
    while(r.hasNext()){
      counter++;
      r.next(event);

      //System.out.println("reading");

      event.read(pt,32100,3);
    }
    return counter;
  }

  


  public static void main(String[] args) {

    GraphErrors gTime = new GraphErrors();
    gTime.attr().setMarkerColor(2);
    gTime.attr().setMarkerSize(16);
    gTime.attr().setTitle("Time vs Number of Threads");
    gTime.attr().setTitleY("Time [s]");
    gTime.attr().setTitleX("Number of Threads");

    GraphErrors gRate = new GraphErrors();
    gRate.attr().setMarkerColor(2);
    gRate.attr().setMarkerSize(16);
    gRate.attr().setTitle("Rate vs Number of Threads");
    gRate.attr().setTitleY("Rate [kHz]");
    gRate.attr().setTitleX("Number of Threads");


    int[] nThreads = new int[]{1,2,4,8}; //,8,12,16};

    double[] times = new double[]{91.333,58.667,46.67,43.667};

    int nEvs=0;

    for (int i=0;i<3;i++){

      int th=nThreads[i];
      

      /* in principle get average and stdev from processor
       * for now just read output out and fill by hand the average time
       */
      double avTime=0;
      /*for(int j=0; j<1;j++){
        Process(th);

      }*/

      //remove later
      avTime=times[i];

      if(nEvs==0){
        nEvs=getNEvs("output_1.h5");
      }

      double avRate=nEvs/avTime;
      avRate=avRate/1000;

      gTime.addPoint(th, avTime, 0, 0);
      gRate.addPoint(th, avRate, 0, 0);

      System.out.printf("\n\n Threads %d, nEvents %d, avTime %.3f, avRate %.3f\n",th,nEvs,avTime,avRate);


    }

    TGCanvas c = new TGCanvas();
    c.setTitle("Time vs nThreads");
    c.draw(gTime);

    TGCanvas c2 = new TGCanvas();
    c2.setTitle("Rate vs nThreads");
    c2.draw(gRate);

  }

}

