
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
public class Level3Converter {
  EJMLModel cf;
  EJMLModel htccer;
  DataList LtoSconv;

  List<float[]> negTracks = new ArrayList<>();
  List<float[]> posTracks = new ArrayList<>();
  List<float[]> all_permutations = new ArrayList<>();

  public Level3Converter() {

  }

  public void load_cf(String path) {
    cf = new EJMLModel(path, ModelType.TANH_LINEAR);
  }

  public void load_htccer(String path) {
    htccer = new EJMLModel(path, ModelType.TANH_LINEAR);
  }


  public void load_LtoStripConv(String path) {
    LtoSconv = DataList.fromCSV("LtoStrip_convTable.csv",
        DataList.range(0, 69), DataList.range(0, 1));
  }

  public Boolean arrayNotEmpty(float[] arr) {
    for (double value : arr) {
      if (value != 0) {
        return true;
      }
    }
    return false;
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

  public int nPart_pSect(List<Level3Candidate> Candidates, int sect) {
    int nPart_pSect = 0;
    for (Level3Candidate part : Candidates) {
      if (part.Cal_Sector == sect) {
        nPart_pSect++;
      }
    }
    return nPart_pSect;
  }

  public int nTrack_pSect(Bank TrackBank, int sect) {
    int nTrack_pSect = 0;
    for (int k = 0; k < TrackBank.getRows(); k++) {
      int pindex = TrackBank.getInt("pindex", k);
      int sectorTrk = TrackBank.getInt("sector", k);
      if (sectorTrk == sect) {
        nTrack_pSect++;
      }
    }
    return nTrack_pSect;
  }

  // 0 1 2 3 4
  // 5 6 7 8 9
  // 10 11 12
  // Bank[] banks = r.getBanks("DC::tdc", "ECAL::adc", "RUN::config", "FTOF::adc",
  // "HTCC::adc",
  // "REC::Particle", "REC::Calorimeter",
  // "REC::Cherenkov","REC::Track","HitBasedTrkg::HBClusters",
  // "HitBasedTrkg::HBTracks","ECAL::clusters", "REC::Scintillator");
  public List<Level3Candidate> getCandidates(Bank[] banks, CompositeNode tr, CompositeNode pt, Boolean isOut) {

    // Instant start_time = Instant.now();

    List<Level3Candidate> ps = new ArrayList<Level3Candidate>();
    for (int i = 0; i < tr.getRows(); i++) {
      // find and initialise Candidates
      Level3Candidate part = new Level3Candidate();
      // part.setShow(true);
      part.find_AI_tracks(tr, pt, i,isOut);
      part.find_ClosestRECParticle_fromPredP(banks[5], 999, 999, 0.2);
      part.read_Cal_Bank(banks[6]);
      part.find_HTCC_ADCs(banks[4]);
      part.read_HTCC_bank(banks[7]);
      part.read_Scint_Bank(banks[12]);
      part.set_pid_label();

      float[] cf_pred = new float[11];
      cf.getOutput(part.track_clusters, cf_pred);
      int[] cf_strips = convLtoStrip(cf_pred);
      part.read_cal_bank_from_cf_pred(cf_pred, cf_strips, banks[1]);

      float[] htcc_pred = new float[1];
      htccer.getOutput(part.get_vars_forhtcc(), htcc_pred);
      part.setPredNphe(htcc_pred[0]*50);

      // part.print();
      // part.printbanks(tr,pt,banks[5],banks[6],banks[4]);

      ps.add(part);
    }

    // Instant finish_time = Instant.now();
    // timeElapsed += Duration.between(start_time, finish_time).toMillis();

    return ps;
  }

  public int hasMatchedPID(List<Level3Candidate> Candidates, int pid) {
    int c = 0;
    for (Level3Candidate part : Candidates) {
      int PID = part.PID;
      if (PID == pid) {
        c++;
      }
    }
    return c;
  }

  public int eventHasRECPID(Bank[] banks, int PID, ArrayList<Level3Particle> parts, CompositeNode tr,
      CompositeNode pt) {
    int c = 0;
    for (int i = 0; i < banks[5].getRows(); i++) {
      Level3Particle part = new Level3Particle();
      part.read_Particle_Bank(i, banks[5]);
      part.find_sector_track(banks[8]);
      part.find_track_clusterIDs(banks[10]);
      part.find_track_clusters(banks[9]);
      part.read_Cal_Bank(banks[6]);
      part.find_sector_cal(banks[6]);
      part.find_HTCC_ADCs(banks[4]);

      if (part.Sector != 0 && part.Track_nSL == 6 && part.Track_chi2 < 350) {
        parts.add(part);

        // part.applyTriangCut();
        if (PID == part.PID) {
          if (PID == 11) { // && arrayNotEmpty(part.HTCC_adcs)
            if (part.Vz < 12 && part.Vz > -13 && Math.abs(part.chi2pid) < 5) { // inbending
            //if(part.Vz<10 && part.Vz>-18 && Math.abs(part.chi2pid)<20){ //outbending
              c++;
            }
          } else {
            c++;
          }
          /*
           * tr.print();
           * pt.print();
           * banks[5].show();
           */

        }

      }
    }
    return c;
  }

  public void convertData(String file, String out, int isPID, int nPart,Boolean isOut) {

    // Specify the path where you want to save the CSV file
    Path filePath = Paths.get(out);

    // Delete the file if it already exists
    try {
      Files.deleteIfExists(filePath);
    } catch (IOException e) {
      e.printStackTrace();
    }

    Boolean notAllFull = true;
    int nEl = 0, nPos = 0, nPip = 0, nPim = 0, nMum = 0, nMup = 0, rEv = 0;
    int nUnnmatchedEl = 0;

    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(filePath))) {
      // String file=dir;

      HipoReader r = new HipoReader(file);

      Event e = new Event();

      // r.getSchemaFactory().show();
      Bank[] banks = r.getBanks("DC::tdc", "ECAL::adc", "RUN::config", "FTOF::adc", "HTCC::adc",
          "REC::Particle", "REC::Calorimeter", "REC::Cherenkov", "REC::Track", "HitBasedTrkg::HBClusters",
          "HitBasedTrkg::HBTracks", "ECAL::clusters", "REC::Scintillator");

      CompositeNode tr = new CompositeNode(3210, 2, "i", 1200);
      CompositeNode pt = new CompositeNode(3210, 3, "i", 1200);

      while (r.hasNext() && notAllFull ) { //&& rEv<10000
        rEv++;

        r.nextEvent(e);
        e.read(banks);

        // System.out.println("DC");
        // banks[0].show();

        /*
         * System.out.println("FTOF");
         * banks[3].show();
         * System.out.println("HTCC");
         * banks[4].show();
         */

        /*
         * System.out.println("REC::Calorimeter");
         * banks[7].show();
         * System.out.println("ECAL::hits");
         * banks[15].show();
         */

        e.read(tr, 32100, 2);
        e.read(pt, 32100, 3);

        List<Level3Candidate> Candidates = getCandidates(banks, tr, pt,isOut);
        // List<Level3Candidate> Candidates = getCandidates_ownc(banks);
        for (Level3Candidate p : Candidates) {
          if (p.unmatched == false && p.Cal_Sector == p.Sector ) { //&& p.Sector==1

            // p.applyTriangCut();
            if (p.PID == 11 && nEl < nPart && p.Charge == -1) {
              if (isPID == 0) {
                if (p.PCALLU_fcf != -2) { // && arrayNotEmpty(p.HTCC_adcs)
                  writer.println(p.get_csv_out_fromCFPred());//_wSF _wHTCCPred_wSF
                  nEl++;
                }
              } else if (isPID == 1) {
                if (!arrayContains(p.CF_out, 0)) {
                  writer.println(p.get_csv_cf_out());
                  nEl++;
                }
              } else if (isPID == 4) {
                if (!arrayContains(p.CF_out_onlyFTOF, 0)) {
                  writer.println(p.get_csv_cf_out_onlyFTOF());
                  nEl++;
                }
              }

            }
            // any neg particle not an e- as BG
            else if (p.PID != 11 && p.Charge == -1 ) {
              if (isPID == 0 && nPim < nPart) {
                if (p.PCALLU_fcf != -2) {
                  writer.println(p.get_csv_out_fromCFPred());
                  nPim++;
                }
              }
            }

            else if (isPID == 2 && p.Charge == 1) {
              if (!arrayContains(p.CF_out, 0)) {
                writer.println(p.get_csv_cf_out());
                nPip++;
              }
            }

            if (isPID == 3 && nEl<nPart && p.Charge == -1) {
              writer.println(p.get_csv_htcc_out());
              nEl++;
            }

          }
        }

        ArrayList<Level3Particle> parts = new ArrayList<Level3Particle>();
        int nRecEl = eventHasRECPID(banks, 11, parts, tr, pt);
        int nMatchedEl = hasMatchedPID(Candidates, 11);
        if (nRecEl > nMatchedEl) {
        
         /*System.out.printf("\nnREC e- %d nMatched e- %d \n",nRecEl,nMatchedEl);
         System.out.println("AI tracks");
         for(int i=0;i<tr.getRows();i++){
         System.out.printf("sector = %2d, charge = %3d segments [ ",
         tr.getInt(1,i),tr.getInt(2,i));
         for(int k = 0; k < 6; k++) System.out.printf(" %9.5f ",tr.getDouble(10+k,i));
         System.out.printf("] p %9.5f  %9.5f %9.5f\n",pt.getDouble(4,i),pt.getDouble(5
         ,i),pt.getDouble(6,i));
         }
         System.out.println("REC tracks");
         for (Level3Particle p:parts){
         if(p.Charge!=0){
         System.out.printf("%d Sector %d charge %d segments [ ",p.PID,p.Sector,
         p.Charge);
         for(int k = 0; k < 6; k++)
         System.out.printf(" %9.5f ",p.track_clusters[k]*112);
         System.out.printf("] p %f %f %f\n",p.Px,p.Py,p.Pz);
         }
         
         }
         banks[5].show();
         banks[8].show();*/
        
          // banks[9].show();

          /*
           * System.out.println("\nall possible tracks");
           * for (float[] ntr:all_permutations){
           * System.out.println(Arrays.toString(ntr));
           * }
           * System.out.println("\nneg pred tracks");
           * for (float[] ntr:negTracks){
           * System.out.println(Arrays.toString(ntr));
           * }
           * System.out.println("pos pred tracks");
           * for (float[] ptr:posTracks){
           * System.out.println(Arrays.toString(ptr));
           * }
           * System.out.println("parts");
           * for (Level3Particle p:parts){
           * System.out.printf("part PID %d Pindex %d Sector %d", p.PID, p.PIndex,
           * p.Sector);
           * System.out.println(Arrays.toString(p.track_clusters));
           * }
           * 
           * System.out.println("REC::Track");
           * banks[8].show();
           * //System.out.println("HB::Cluster");
           * //banks[9].show();
           * System.out.println("HB::Tracks");
           * banks[10].show();
           */

          nUnnmatchedEl += nRecEl - nMatchedEl;
        }

        if (nEl >= nPart && nPim >= nPart) { // && nPip>=nPart ){// && nPos>= nPart && nMum>=nPart && nMup>=
                                             // nPart){
          notAllFull = false;
        }

        if(isPID==1 && nEl >= nPart){
          notAllFull = false;
        }

        if(isPID==2 && nPip >= nPart){
          notAllFull = false;
        }

      }

    } catch (IOException e) {
      e.printStackTrace();
    }

    double eff = ((float) nEl) / (((float) nEl) + ((float) nUnnmatchedEl));
    System.out.printf("Found %d 11, %d -11, %d -211, %d 211, %d -13, %d 13\n", nEl, nPos, nPim, nPip, nMum, nMup);
    System.out.printf("Nb UnMatched 11: %d, eff %f\n", nUnnmatchedEl, eff);
    System.out.printf("Read %d events\n\n", rEv);

  }

  public static void main(String[] args) {

    // to run, in order
    // /open Level3Particle.java
    // /open Level3Candidate.java
    // /open Level3Converter.java
    // Level3Converter.main(new String[]{});

    String file = "/Users/tyson/data_repo/trigger_data/rgd/018326/run_18326_1_wAIBanks.h5";
    String filet = "/Users/tyson/data_repo/trigger_data/rgd/018326/run_18326_3_wAIBanks.h5";
    String outDir = "/Users/tyson/data_repo/trigger_data/rgd/018326/for_caos_pid/";
    Boolean isOut=false;

    //String file = "/Users/tyson/data_repo/trigger_data/rgd/018777/run_18777_1_wAIBanks.h5";
    //String filet = "/Users/tyson/data_repo/trigger_data/rgd/018777/run_18777_3_wAIBanks.h5";
    //String outDir = "/Users/tyson/data_repo/trigger_data/rgd/018777/for_caos_pid/";
    //Boolean isOut=true;

    String field="";
    if(isOut){field="_outbending";}

    Level3Converter conv = new Level3Converter();
    conv.load_cf("cf_el_wFTOF"+field+".network");
    conv.load_LtoStripConv("LtoStrip_convTable.csv");
    conv.load_htccer("htccer_allNeg"+field+".network");

    // isPID 0 - training for pid, 1 - training for negative particle cf, 
    //2 - training for pos part cf, 3 - training HTCC, 4 - only FTOF

    //conv.convertData(file,outDir+"train_fromcfpred_allNegBG.csv",0,100000,isOut);//50000 //_wHTCCPred_wSF
    //conv.convertData(filet,outDir+"test_fromcfpred_allNegBG.csv",0,50000,isOut);//50000 //_wHTCCPred_wSF

    //conv.convertData(file, outDir + "train_cf_wFTOF.csv", 1, 100000,isOut);// 50000
    //conv.convertData(filet, outDir + "test_cf_wFTOF.csv", 1, 50000,isOut);// 50000

    conv.convertData(file, outDir + "train_cf_onlyFTOF.csv", 4, 100000,isOut);// 50000
    conv.convertData(filet, outDir + "test_cf_onlyFTOF.csv", 4, 50000,isOut);// 50000

    //conv.convertData(file,outDir+"train_allNeg_htcc.csv",3,100000,isOut);//50000
    //conv.convertData(filet,outDir+"test_allNeg_htcc.csv",3,10000,isOut);//50000

  }
}
