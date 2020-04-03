class AminoAcidLL{
  char aminoAcid;
  String[] codons;
  int[] counts;
  AminoAcidLL next;

  AminoAcidLL(){

  }


  /********************************************************************************************/
  /* Creates a new node, with a given amino acid/codon 
   * pair and increments the codon counter for that codon.
   * NOTE: Does not check for repeats!! */
  AminoAcidLL(String inCodon){

    this.aminoAcid = AminoAcidResources.getAminoAcidFromCodon(inCodon);
    this.codons = AminoAcidResources.getCodonListForAminoAcid(this.aminoAcid);
    this.counts = new int[codons.length];
    this.next = null;

  }

  /********************************************************************************************/
  /* Recursive method that increments the count for a specific codon:
   * If it should be at this node, increments it and stops, 
   * if not passes the task to the next node. 
   * If there is no next node, add a new node to the list that would contain the codon. 
   */
  private void addCodon(String inCodon){

    if(aminoAcid == AminoAcidResources.getAminoAcidFromCodon(inCodon)) {
      for (int i = 0; i < codons.length; i++) {
        if (codons[i].equals(inCodon))
          counts[i]++;
      }
    }else if(next != null){
      next.addCodon(inCodon);
    }else{
      next = new AminoAcidLL(inCodon);
      for (int i = 0; i < codons.length; i++) {
        if (codons[i].equals(inCodon))
          counts[i]++;
      }
    }

  }

  /********************************************************************************************/
  /* Shortcut to find the total number of instances of this amino acid */
  private int totalCount(){

    int sum = 0;
    for(int i = 0; i < counts.length; i++)
        sum += counts[i];
    return sum;

  }

  /********************************************************************************************/
  /* helper method for finding the list difference on two matching nodes
  *  must be matching, but this is not tracked */
  private int totalDiff(AminoAcidLL inList){

    return Math.abs(totalCount() - inList.totalCount());

  }


  /********************************************************************************************/
  /* helper method for finding the list difference on two matching nodes
  *  must be matching, but this is not tracked */
  private int codonDiff(AminoAcidLL inList){
    int diff = 0;
    for(int i=0; i < codons.length; i++){
      diff += Math.abs(counts[i] - inList.counts[i]);
    }
    return diff;
  }

  /********************************************************************************************/
  /* Recursive method that finds the differences in **Amino Acid** counts. 
   * the list *must* be sorted to use this method */
  public int aminoAcidCompare(AminoAcidLL inList){

    if(!inList.isSorted())
      inList.sort(inList);

    if(!this.isSorted())
      this.sort(this);

    return aminoHelper(inList, 0);

  }

  public int aminoHelper(AminoAcidLL inList, int diff){

    if(this == null && inList == null)
      return diff;

    while(this == null || aminoAcid > inList.aminoAcid){
      diff += inList.totalCount();
      inList = inList.next;
    }

    while(inList == null || aminoAcid < inList.aminoAcid){
      diff += totalCount();
      diff += next.aminoHelper(inList, diff);
    }

    if(aminoAcid == inList.aminoAcid){
      diff += totalDiff(inList);
      diff += next.aminoHelper(inList.next, diff);
    }
    return diff;
  }

  /********************************************************************************************/
  /* Same ad above, but counts the codon usage differences
   * Must be sorted. */
  public int codonCompare(AminoAcidLL inList){

    if(!inList.isSorted())
      inList.sort(inList);

    if(!this.isSorted())
      this.sort(this);

    return codonHelper(inList, 0);

  }

  public int codonHelper(AminoAcidLL inList, int diff){

    if(this == null && inList == null)
      return diff;

    while(this == null || aminoAcid > inList.aminoAcid){
      diff += inList.totalCount();
      inList = inList.next;
    }

    while(inList == null || aminoAcid < inList.aminoAcid){
      diff += totalCount();
      diff += next.codonHelper(inList, diff);
    }

    if(aminoAcid == inList.aminoAcid){
      diff += codonDiff(inList);
      diff += next.codonHelper(inList.next, diff);
    }
    return diff;
  }


  /********************************************************************************************/
  /* Recursively returns the total list of amino acids in the order that they are in in the linked list. */
  public char[] aminoAcidList(){

    char[] letters = new char[size()];
    aminoAcidListHelper(letters, 0);
    return letters;

  }

  public void aminoAcidListHelper(char[] letters, int index){

    letters[index] = aminoAcid;
    if(next == null)
      return;
    next.aminoAcidListHelper(letters, index + 1);

  }

  public int size(){

    if(next == null)
      return 1;
    return 1 + next.size();

  }

  /********************************************************************************************/
  /* Recursively returns the total counts of amino acids in the order that they are in in the linked list. */
  public int[] aminoAcidCounts(){

    int[] totalAmino = new int[size()];
    aminoAcidCountsHelper(totalAmino, 0);
    return totalAmino;

  }

  public void aminoAcidCountsHelper(int[] totalAmino, int index){

    totalAmino[index] = totalCount();
    if(next == null)
      return;
    next.aminoAcidCountsHelper(totalAmino, index + 1);

  }



  /********************************************************************************************/
  /* recursively determines if a linked list is sorted or not */
  public boolean isSorted(){

    if(next == null)
      return true;
    if(aminoAcid > next.aminoAcid)
      return false;
    if(next.next != null)
      return next.isSorted();

    return true;

  }


  /********************************************************************************************/
  /* Static method for generating a linked list from an RNA sequence */
  public static AminoAcidLL createFromRNASequence(String inSequence) {

    if(inSequence.length() < 3)
      return null;

    String inCodon = inSequence.substring(0, 3);

    if(inCodon.equals("UGA") || inCodon.equals("UAG") || inCodon.equals("UAA"))
      return null;

    AminoAcidLL inList = new AminoAcidLL(inCodon);
    return createHelper(inSequence, inList);

  }

  public static AminoAcidLL createHelper(String inSequence, AminoAcidLL inListStart){

    AminoAcidLL inList = inListStart;

    if(inSequence.length() < 3)
      return null;

    String inCodon = inSequence.substring(0, 3);
    char aminoAcid = AminoAcidResources.getAminoAcidFromCodon(inCodon);

    if(inCodon.equals("UGA") || inCodon.equals("UAG") || inCodon.equals("UAA"))
      return null;

    inList.addCodon(inCodon);

    return createHelper(inSequence.substring(3), inListStart);

  }



  /********************************************************************************************/
  /* sorts a list by amino acid character*/
  public static AminoAcidLL sort(AminoAcidLL inList){

    AminoAcidLL inListStart = inList;

    AminoAcidLL resultHead = new AminoAcidLL();
    resultHead.next = inList;

    while(inList.next != null) {

      inList = inList.next;
      AminoAcidLL previous = resultHead;
      AminoAcidLL current = resultHead.next;

      while (inList.aminoAcid > current.aminoAcid && current != null) {
        previous = previous.next;
        current = current.next;

      }
      previous.next = inList;

    }
    return resultHead.next;
  }
}