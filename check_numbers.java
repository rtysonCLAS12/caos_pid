HipoReader r = new HipoReader("output_test_inbending_large_newNetwork.h5");

System.out.println("read");


CompositeNode tr = new CompositeNode(32100,2,"i",1200);
CompositeNode pt = new CompositeNode(32100,3,"i",1200);
CompositeNode cf = new CompositeNode(32100,11,"i",1200);
CompositeNode pt_test = new CompositeNode(32100,22,"i",1200);
CompositeNode cf_test = new CompositeNode(32100,21,"i",24);
CompositeNode pt_test_nAI = new CompositeNode(32100,32,"i",1200);
CompositeNode cf_test_nAI = new CompositeNode(32100,31,"i",24);

Event event = new Event();
int counter = 0;
int hasInstaEl_hasAITrack=0,hasInstaEl_hasAIEl=0,hasInstaEl_hasTrack=0,hasInstaEl_hasEl=0, hasInstaEl=0;
int hasAIEl_hasInstaTrack,hasAIEl_hasInstaEl=0,hasAIEl_hasEl=0,hasAIEl_hasTrack=0,hasAIEl=0;
int hasEl_hasInstaTrack=0,hasEl_hasInstaEl=0,hasEl_hasAIEl=0,hasEl_hasAITrack=0,hasEl=0;
while(r.hasNext()){
    counter++;
    r.next(event);

    //System.out.println("reading");
    event.read(tr,32100,2);
    //System.out.println("reading 2");
    event.read(pt,32100,3);
    //System.out.println("reading 3");
    event.read(cf,32100,11);
    //System.out.println("reading 4");
    event.read(pt_test,32100,22);
    //System.out.println("reading 5");
    event.read(cf_test,32100,21);
    //System.out.println("reading 6");
    event.read(pt_test_nAI,32100,32);
    //System.out.println("reading 7");
    event.read(cf_test_nAI,32100,31);
    //System.out.println("read");
    int nrows = tr.getRows();
    int nrowstest=pt_test.getRows();
    int nrowstest_nAI=pt_test_nAI.getRows();

    //check if InstaRec has e-
    for(int j = 0; j < nrows; j++){
      if(tr.getInt(1,j)==1){
        //has Insta El
        if(pt.getInt(3,j) == 11){
          hasInstaEl++;
          double Pz_o=pt.getDouble(6,j);

          for(int i = 0; i < nrowstest; i++){
            if(pt_test.getInt(2,i)==1){
              double resPz=Math.abs(pt_test.getDouble(6,i)-Pz_o)/Math.abs(Pz_o);
              if(pt.getInt(1,j)==pt_test.getInt(1,i) && resPz<0.2){
                hasInstaEl_hasAITrack++;
              }
              if(pt_test.getInt(3,i)==11){
                hasInstaEl_hasAIEl++;
              }
            }
          }

          for(int i = 0; i < nrowstest_nAI; i++){
            if(pt_test_nAI.getInt(2,i)==1){
              double resPz=Math.abs(pt_test_nAI.getDouble(6,i)-Pz_o)/Math.abs(Pz_o);
              if(pt.getInt(1,j)==pt_test_nAI.getInt(1,i) && resPz<0.2){
                hasInstaEl_hasTrack++;
              }
              if(pt_test_nAI.getInt(3,i)==11){
                hasInstaEl_hasEl++;
              }
            }
          }
          
        }
      }
    }

    //has REC AI e-
    for(int j = 0; j < nrowstest; j++){
      if(pt_test.getInt(2,j)==1 ){
        if(pt_test.getInt(3,j)==11 ){
          hasAIEl++;

          double Pz_o=pt_test.getDouble(6,j);

          for(int i = 0; i < nrows; i++){
            if(tr.getInt(1,i)==1){
              double resPz=Math.abs(pt.getDouble(6,i)-Pz_o)/Math.abs(Pz_o);
              if(pt.getInt(1,i)==pt_test.getInt(1,j) && resPz<0.2){
                hasAIEl_hasInstaTrack++;
              }
              if(pt.getInt(3,i) == 11){
                hasAIEl_hasInstaEl++;
              }
            }
          }

          for(int i = 0; i < nrowstest_nAI; i++){
            if(pt_test_nAI.getInt(2,i)==1){
              double resPz=Math.abs(pt_test_nAI.getDouble(6,i)-Pz_o)/Math.abs(Pz_o);
              if(pt_test.getInt(1,j)==pt_test_nAI.getInt(1,i) && resPz<0.2){
                hasAIEl_hasTrack++;
              }
              if(pt_test_nAI.getInt(3,i)==11){
                hasAIEl_hasEl++;
              }
            }
          }
          
        }
      }
    }

    //has AI e-
    for(int j = 0; j < nrowstest_nAI; j++){
      if(pt_test_nAI.getInt(2,j)==1){
        if(pt_test_nAI.getInt(3,j)==11){
          hasEl++;

          double Pz_o=pt_test_nAI.getDouble(6,j);

          for(int i = 0; i < nrows; i++){
            if(tr.getInt(1,i)==1){
              double resPz=Math.abs(pt.getDouble(6,i)-Pz_o)/Math.abs(Pz_o);
              if(pt.getInt(1,i)==pt_test_nAI.getInt(1,j) && resPz<0.2){
                hasEl_hasInstaTrack++;
              }
              if(pt.getInt(3,i) == 11){
                hasEl_hasInstaEl++;
              }
            }
          }

          for(int i = 0; i < nrowstest; i++){
            if(pt_test.getInt(2,i)==1){
              double resPz=Math.abs(pt_test.getDouble(6,i)-Pz_o)/Math.abs(Pz_o);
              if(pt_test.getInt(1,i)==pt_test_nAI.getInt(1,j) && resPz<0.2){
                hasEl_hasAITrack++;
              }
              if(pt_test_nAI.getInt(3,i)==11){
                hasEl_hasAIEl++;
              }
            }
          }
          
        }
      }
    }
}

System.out.printf("\n InstaRec has e-: \t Total %d \n",hasInstaEl);
System.out.printf("\t REC AI has e- %d \n",hasInstaEl_hasAIEl);
System.out.printf("\t REC AI has same track %d \n",hasInstaEl_hasAITrack);
System.out.printf("\t REC has e- %d \n",hasInstaEl_hasEl);
System.out.printf("\t REC has same track %d \n",hasInstaEl_hasTrack);

System.out.printf("\n\n REC AI has e-:\t Total %d \n",hasAIEl);
System.out.printf("\t InstaRec has e- %d \n",hasAIEl_hasInstaEl);
System.out.printf("\t InstaRec has same track %d \n",hasAIEl_hasInstaTrack);
System.out.printf("\t REC has e- %d \n",hasAIEl_hasEl);
System.out.printf("\t REC has same track %d \n",hasAIEl_hasTrack);

System.out.printf("\n\n REC has e-:\t Total %d \n",hasEl);
System.out.printf("\t InstaRec has e- %d \n",hasEl_hasInstaEl);
System.out.printf("\t InstaRec has same track %d \n",hasEl_hasInstaTrack);
System.out.printf("\t REC AI has e- %d \n",hasEl_hasAIEl);
System.out.printf("\t REC AI has same track %d \n",hasEl_hasAITrack);

/exit
