HipoReader r = new HipoReader("output_test_inbending_large.h5");

System.out.println("read");


CompositeNode pt = new CompositeNode(32100,3,"i",1200);

Event event = new Event();
int counter = 0;
double nAllEls=0,nElAndTrack=0,nTot=0;

while(r.hasNext()){
  counter++;
  r.next(event);

  event.read(pt,32100,3);
  int nrows = pt.getRows();

  List<Integer> elSectors = new ArrayList<Integer>();
  List<Integer> notElSectors = new ArrayList<Integer>();  

  for(int j = 0; j < nrows; j++){
    int charge=pt.getInt(1,j);
    int sector=pt.getInt(2,j);
    int pid=pt.getInt(3,j);

    //System.out.printf("bla %d %d %d\n",charge,sector,pid);

    if(charge==-1){
      if(pid==11){
        elSectors.add(sector);
      } else{
        notElSectors.add(sector);
      }
    }
  }

  if(elSectors.size()>0){
    for(int i=0;i<elSectors.size();i++){
      for (int k = i + 1; k < elSectors.size(); k++) {
        if(elSectors.get(i)==elSectors.get(k)){
          nAllEls++;
          nTot++;
        }
      }
  
      if(notElSectors.size()>0){
        for (int k = 0; k < notElSectors.size(); k++) {
          if(elSectors.get(i)==notElSectors.get(k)){
            nElAndTrack++;
            nTot++;
          }
        }
      }
    }
  }

}

double r1=nAllEls/nTot;
double r2=nElAndTrack/nTot;

System.out.printf("Have %f events with two or more neg tracks in same sector\n",nTot);
System.out.printf("Of these, %f have two electrons eg ratio %f \n",nAllEls,r1);
System.out.printf("Of these, %f have electron and other pid eg ratio %f \n",nElAndTrack,r2);

/exit
