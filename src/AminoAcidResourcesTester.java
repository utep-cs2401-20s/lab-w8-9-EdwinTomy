import static org.junit.jupiter.api.Assertions.*;
import java.util.concurrent.TimeUnit;

import org.junit.jupiter.api.Test;

class AminoAcidResourcesTester{

  @Test
  public void allCodons(){
    char[] rna = {'A','C','U','G'};
    char[] aa = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W'};
    for(int i=0; i<4; i++){
      for(int j=0; j<4; j++){
        for(int k=0; k<4;k++){
          String s = new String(new char[]{rna[i],rna[j],rna[k]});
          char aaOut = AminoAcidResources.getAminoAcidFromCodon(s);
          if(aaOut != '*'){
            String[] codonList = AminoAcidResources.getCodonListForAminoAcid(aaOut);
            boolean found = false;
            for(int l=0; l<codonList.length; l++){
              found |= (codonList[l].equals(s));
            }
            if(!found) System.err.println("Codon " + s + " not found, said AA was " + aaOut);
          }

          aaOut = AminoAcidResources.getAminoAcidFromCodon(s.toLowerCase());
          if(aaOut != '*'){
            String[] codonList = AminoAcidResources.getCodonListForAminoAcid(aaOut);
            boolean found = false;
            for(int l=0; l<codonList.length; l++){
              found |= (codonList[l].equals(s));
            }
            if(!found) System.err.println("Codon " + s + " not found, said AA was " + aaOut);
          }
        }
      }
    }
  }

  @Test
  public void allAAs(){

    char[] aa = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W'};
    for(int i=0; i<aa.length; i++){
      String[] codonList = AminoAcidResources.getCodonListForAminoAcid(aa[i]);
      for(int l=0; l<codonList.length; l++){
        if(aa[i] != AminoAcidResources.getAminoAcidFromCodon(codonList[l])){
          System.err.println("AA " + aa[i] + " not found, said codon was " + codonList[l]);
        }
      }

      codonList = AminoAcidResources.getCodonListForAminoAcid(Character.toLowerCase(aa[i]));
      for(int l=0; l<codonList.length; l++){
        if(aa[i] != AminoAcidResources.getAminoAcidFromCodon(codonList[l])){
          System.err.println("AA " + aa[i] + " not found, said codon was " + codonList[l]);
        }
      }
    }
  }

  //Test case for constructor: PASSED
  @Test
  public void constructorTest(){

    System.out.println("Test case for AminoAcidLL()");
    AminoAcidLL amino = new AminoAcidLL("GCC");
    String[] expected = {"GCG","GCA","GCC","GCU"};

    assertArrayEquals(expected, amino.codons);

  }

  //The following test cases also test for createFromRNASequence();
  //Test case for aminoAcidList(): PASSED
  @Test
  public void aminoAcidListTest(){

    System.out.println("Test case for aminoAcidList()");
    AminoAcidLL head = AminoAcidLL.createFromRNASequence("CCGAUGAAC");


    assertArrayEquals(new char[]{'P', 'M', 'N'}, head.aminoAcidList());

  }

  //Test case for aminoAcidsCounts(): PASSED
  @Test
  public void aminoAcidCountsTest(){

    System.out.println("Test case for aminoAcidCounts()");
    AminoAcidLL head = AminoAcidLL.createFromRNASequence("CCGAUGAACCCAAAUCCCCUGCUA");

    assertArrayEquals(new int[]{3, 1, 2, 2}, head.aminoAcidCounts());

  }

  //Test case for aminoAcidCompare(): PASSED
  @Test
  public void aminoAcidCompareTest(){

    System.out.println("Test case for aminoAcidCompare()");
    AminoAcidLL headOne = AminoAcidLL.createFromRNASequence("UGCGACGAGUUC");
    AminoAcidLL headTwo = AminoAcidLL.createFromRNASequence("GCAUGUGAUGAAUUUGGA");

    assertEquals(2, headOne.aminoAcidCompare(headTwo));

  }

  //Test case for codonCompare(): PASSED
  @Test
  public void codonCompareTest(){

    System.out.println("Test case for codonCompare()");
    AminoAcidLL headOne = AminoAcidLL.createFromRNASequence("UGCGACGAGUUC");
    AminoAcidLL headTwo = AminoAcidLL.createFromRNASequence("GCAUGUGAUGAAUUUGGA");

    assertEquals(10, headOne.codonCompare(headTwo));

  }

  //The following test cases also test for createFromRNASequence() with stop condons;
  //Test case for aminoAcidList(): PASSED
  @Test
  public void aminoAcidListTestX(){

    System.out.println("Test case for aminoAcidList() with stop");
    AminoAcidLL head = AminoAcidLL.createFromRNASequence("CCGAUGAACUGACCG");


    assertArrayEquals(new char[]{'P', 'M', 'N'}, head.aminoAcidList());

  }

  //Test case for aminoAcidsCounts(): PASSED
  @Test
  public void aminoAcidCountsTestX(){

    System.out.println("Test case for aminoAcidCounts() with stop");
    AminoAcidLL head = AminoAcidLL.createFromRNASequence("CCGAUGAACCCAAAUCCCCUGCUAUAGCCGCCGCCG");

    assertArrayEquals(new int[]{3, 1, 2, 2}, head.aminoAcidCounts());

  }

  //Test case for aminoAcidCompare(): PASSED
  @Test
  public void aminoAcidCompareTestX(){

    System.out.println("Test case for aminoAcidCompare() with stop");
    AminoAcidLL headOne = AminoAcidLL.createFromRNASequence("UGCGACGAGUUCUAAUGCUGC");
    AminoAcidLL headTwo = AminoAcidLL.createFromRNASequence("GCAUGUGAUGAAUUUGGAUAAGGAGGA");

    assertEquals(2, headOne.aminoAcidCompare(headTwo));

  }

  //Test case for codonCompare(): PASSED
  @Test
  public void codonCompareTestX(){

    System.out.println("Test case for codonCompare() with stop");
    AminoAcidLL headOne = AminoAcidLL.createFromRNASequence("UGCGACGAGUUCUAAUGC");
    AminoAcidLL headTwo = AminoAcidLL.createFromRNASequence("GCAUGUGAUGAAUUUGGAUAAUGU");

    assertEquals(10, headOne.codonCompare(headTwo));

  }

  //Test case for isSorted(): PASSED
  @Test
  public void isSortedTest(){

    System.out.println("Test case for isSorted()");
    AminoAcidLL first = new AminoAcidLL("GCG");
    AminoAcidLL second = new AminoAcidLL("GAC");
    AminoAcidLL third = new AminoAcidLL("UGC");

    first.next = second;
    second.next = third;

    assertEquals(false, first.isSorted());

  }

  //Test case for sort(): PASSED
  @Test
  public void sortTest(){

    System.out.println("Test case for sort()");
    AminoAcidLL first = new AminoAcidLL("GCG");
    AminoAcidLL second = new AminoAcidLL("GAC");
    AminoAcidLL third = new AminoAcidLL("UGC");

    first.next = second;
    second.next = third;
    first.sort(first);

    assertEquals(true, first.isSorted());

  }

}