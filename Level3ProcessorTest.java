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

/**
 *
 * @author tyson
 */
public class Level3ProcessorTest {

  public Level3ProcessorTest() {

  }

  public static void main(String[] args) {
    //String file = "/Users/tyson/data_repo/trigger_data/rgd/018326/run_18326_1_wAIBanks.h5";
    String file = "/Users/tyson/data_repo/trigger_data/rgd/018777/run_18777_3_wAIBanks.h5";


    DataStream<HipoReader, HipoWriter, Event> str = new DataStream<HipoReader, HipoWriter, Event>();
    str.show();

    DataFrame<Event> frame = new DataFrame<>();
    HipoReader source = new HipoReader();
    HipoWriter sync = new HipoWriter();

    sync.open("output_RecEl_noInstaEl.h5");

    source.open(file);

    DataWorker<HipoReader, Event> worker = new DataWorker<HipoReader, Event>() {

      public final AtomicInteger counter = new AtomicInteger(0);
      Level3PIDUtils utils = new Level3PIDUtils();
      private Schema ECALadc_schema = null;
      private Schema HTCCadc_schema = null;
      private Schema ECALrec_schema = null;
      private Schema SCINTrec_schema = null;
      private Schema PARTrec_schema = null;
      private Schema TRACKrec_schema = null;
      EJMLModel cf;
      DataList LtoSconv;
      EJMLModel pider;

      @Override
      public void execute(Event e) {

        if (counter.intValue() < 10000) {

          // input HTCC ECAL ADC
          Bank bHTCC = new Bank(HTCCadc_schema);
          Bank bECAL = new Bank(ECALadc_schema);
          Bank bSCINTrec = new Bank(SCINTrec_schema);
          Bank bECALrec = new Bank(ECALrec_schema);
          Bank bPARTrec = new Bank(PARTrec_schema);
          Bank bTRACKrec = new Bank(TRACKrec_schema);
          e.read(bHTCC);
          e.read(bECAL);
          e.read(bSCINTrec);
          e.read(bECALrec);
          e.read(bPARTrec);
          e.read(bTRACKrec);

          // input track
          CompositeNode tr = new CompositeNode(32100, 2, "i", 1200);
          CompositeNode pt = new CompositeNode(32100, 3, "i", 1200);
          e.read(tr, 32100, 2);
          e.read(pt, 32100, 3);

          // output
          CompositeNode bcf_out = new CompositeNode(32100, 11, "sf9f", 24);
          CompositeNode bcf_out_test = new CompositeNode(32100, 21, "sf9f", 0);
          CompositeNode pt_out = new CompositeNode(32100, 3, "sssifffffffs", 1200);
          CompositeNode pt_out_test = new CompositeNode(32100, 22, "sssifffffffs", 0);

          int nInstaEl=0;
          for (int i = 0; i < tr.getRows(); i++) {

            int Sector = tr.getInt(1, i);
            int Charge = -tr.getInt(2, i); // !!!!!!!!!!!!!!!!!!! for outbending data flip charge
            float[] track_clusters = new float[6];

            for (int k = 0; k < 6; k++) {
              track_clusters[k] = (float) tr.getDouble(10 + k, i) / 112;
            }

            float[] cf_pred = new float[11];
            cf.getOutput(track_clusters, cf_pred);
            int[] cf_strips = new int[] { -2, -2, -2, -2, -2, -2, -2, -2, -2 };
            utils.convLtoStrip(cf_pred, cf_strips, LtoSconv);

            int outpid = 211;
            float resp = 0;
            if (Charge == -1) {
              float[] vars_pid = new float[35];
              utils.get_varsPID(cf_pred, cf_strips, track_clusters, vars_pid, bECAL, bHTCC, Sector);
              float[] pid_pred = new float[2];
              pider.getOutput(vars_pid, pid_pred);
              resp = pid_pred[0];
              if (pid_pred[0] > 0.1) {
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

            if(Sector==1 && outpid==11){
              nInstaEl++;
            }

          }

          int nrows = 0, nRECel=0;
          for (int i = 0; i < bPARTrec.getRows(); i++) {

            int[] cf_out_test = new int[10];
            float[] apt_out_test = new float[12];
            double track_chi2 = utils.readTestBanks(i, bPARTrec, bECALrec, bSCINTrec, bTRACKrec, cf_out_test,
                apt_out_test, LtoSconv);

            if (Math.abs(apt_out_test[11]) >= 2000 && Math.abs(apt_out_test[11]) < 4000 && track_chi2<350) {
              bcf_out_test.setRows(nrows + 1);
              bcf_out_test.putShort(0, nrows, (short) i);
              // Predicted ECAL Position
              for (int j = 0; j < 9; j++) {
                bcf_out_test.putFloat(j + 1, nrows, cf_out_test[j]);
              }
              // Predicted FTOF Position
              bcf_out_test.putFloat(10, nrows, cf_out_test[9]);

              pt_out_test.setRows(nrows + 1);
              for (int j = 0; j < 3; j++) {
                pt_out_test.putShort(j, nrows, (short) apt_out_test[j]);
              }
              pt_out_test.putInt(3, nrows, (int) apt_out_test[3]);
              // rest of bank is the same
              for (int j = 4; j < 11; j++) {
                pt_out_test.putFloat(j, nrows, apt_out_test[j]);
              }
              pt_out_test.putShort(11, nrows, (short) apt_out_test[11]);
              nrows++;

              if(apt_out_test[3]==11 && apt_out_test[2]==1){
                nRECel++;
              }

            }
          }

          // below fails with:
          // error replacing node ( 32100, 3) dues to size difference 126,132
          // e.replace(32100, 3, pt_out);

          // Messing with outstream because e.remove prints LENGTH = somelength
          // this was annoying me

          //if(nRECel>0 && nInstaEl==0){
            PrintStream original = System.out;
            System.setOut(new PrintStream(OutputStream.nullOutputStream()));
  
            e.remove(32100, 3);
  
            System.setOut(original);
            e.write(pt_out);
  
            e.write(bcf_out); // 32100,11
  
            e.write(bcf_out_test);
            e.write(pt_out_test);
  
            counter.set(counter.intValue() + 1);

          //}
         

          try {
            Thread.sleep(25);
          } catch (InterruptedException ex) {
            Logger.getLogger(DataStream.class.getName()).log(Level.SEVERE, null, ex);
          }
        }

      }

      @Override
      public boolean init(HipoReader src) {
        // ECAL ADC schema used for prediction
        ECALadc_schema = src.getSchemaFactory().getSchema("ECAL::adc");
        HTCCadc_schema = src.getSchemaFactory().getSchema("HTCC::adc");
        ECALrec_schema = src.getSchemaFactory().getSchema("REC::Calorimeter");
        SCINTrec_schema = src.getSchemaFactory().getSchema("REC::Scintillator");
        PARTrec_schema = src.getSchemaFactory().getSchema("REC::Particle");
        TRACKrec_schema = src.getSchemaFactory().getSchema("REC::Track");

        // load networks
        cf = new EJMLModel("cf_el_wFTOF_outbending.network", ModelType.TANH_LINEAR);
        LtoSconv = DataList.fromCSV("LtoStrip_convTable.csv",
            DataList.range(0, 69), DataList.range(0, 1));
        pider = new EJMLModel("pid_elNegBG_fromcfpred_outbending.network", ModelType.SOFTMAX);
        return true;
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
    str.threads(4);
    str.withSource(source).withOutput(sync).withFrame(frame).consumer(worker).run();

    str.show();
  }

}