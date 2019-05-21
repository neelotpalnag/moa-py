//  SolutionSet.Java
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>
//
//  Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
// 
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package jmetal.core;

import java.io.*;
import java.sql.Time;
import java.util.*; 

import jmetal.init.Params;
import jmetal.qualityIndicator.GeneralizedSpread;
import jmetal.qualityIndicator.Hypervolume;
import jmetal.qualityIndicator.Spread;
import jmetal.util.Configuration;
import jmetal.util.JMException;
import jmetal.util.Ranking;

/** 
 * Class representing a SolutionSet (a set of solutions)
 */
public class SolutionSet implements Serializable {
    
  /**
   * Stores a list of <code>solution</code> objects.
   */
  protected List<Solution> solutionsList_;
  
  /** 
   * Maximum size of the solution set 
   */
  private int capacity_ = 0; 
    
  /**
   * Constructor.
   * Creates an unbounded solution set.
   */
  public SolutionSet() {
    solutionsList_ = new ArrayList<Solution>();
  } // SolutionSet
    
  /** 
   * Creates a empty solutionSet with a maximum capacity.
   * @param maximumSize Maximum size.
   */
  public SolutionSet(int maximumSize){    
    solutionsList_ = new ArrayList<Solution>();
    capacity_      = maximumSize;
  } // SolutionSet
  
  /**
   * Inserts a new solution into the SolutionSet
   * @param solution The <code>Solution</code> to store
   */
  public void addToAll(Solution solution){
      solutionsList_.add(solution);
  }
  
 /** 
  * Inserts a new solution into the SolutionSet. 
  * @param solution The <code>Solution</code> to store
  * @return True If the <code>Solution</code> has been inserted, false 
  * otherwise. 
  */
  public boolean add(Solution solution) {
    if (solutionsList_.size() == capacity_) {
      Configuration.logger_.severe("The population is full");
      Configuration.logger_.severe("Capacity is : "+capacity_);
      Configuration.logger_.severe("\t Size is: "+ this.size());
      return false;
    } // if
    
    solutionsList_.add(solution);
    return true;
  } // add
    
  /**
   * Returns the ith solution in the set.
   * @param i Position of the solution to obtain.
   * @return The <code>Solution</code> at the position i.
   * @throws IndexOutOfBoundsException.
   */
  public Solution get(int i) {
    if (i >= solutionsList_.size()) {
      throw new IndexOutOfBoundsException("Index out of Bound "+i);
    }
    return solutionsList_.get(i);
  } // get
      
  /**
   * Returns the maximum capacity of the solution set
   * @return The maximum capacity of the solution set
   */
  public int getMaxSize(){
    return capacity_ ;
  } // getMaxSize
    
  
  /**
   * by Miha
   * Returns the maximum capacity of the solution set
   * @return The maximum capacity of the solution set
   */
  public void setMaxSize(int newCapacity){
    capacity_  = newCapacity;
  } // setMaxSize
  
  
  /** 
   * Sorts a SolutionSet using a <code>Comparator</code>.
   * @param comparator <code>Comparator</code> used to sort.
   */
  public void sort(Comparator comparator){
    if (comparator == null) {
      Configuration.logger_.severe("No criterium for compare exist");
      return ;
    } // if
    Collections.sort(solutionsList_,comparator);
  } // sort
      

  /** 
   * Returns the index of the best Solution using a <code>Comparator</code>.
   * If there are more than one occurrences, only the index of the first one is returned
   * @param comparator <code>Comparator</code> used to compare solutions.
   * @return The index of the best Solution attending to the comparator or 
   * <code>-1<code> if the SolutionSet is empty
   */
  public int indexBest(Comparator comparator){
    
   
    if ((solutionsList_ == null) || (this.solutionsList_.isEmpty())) {
        return -1;
    }
    
    int index = 0; 
    Solution bestKnown = solutionsList_.get(0), candidateSolution;
    int flag;
    for (int i = 1; i < solutionsList_.size(); i++) {        
        candidateSolution = solutionsList_.get(i);
        flag = comparator.compare(bestKnown, candidateSolution);
        if (flag == -1) {
            index = i;
            bestKnown = candidateSolution; 
        }
    }
    
    return index;
        
  } // indexBest
  
  
  /** 
   * Returns the best Solution using a <code>Comparator</code>.
   * If there are more than one occurrences, only the first one is returned
   * @param comparator <code>Comparator</code> used to compare solutions.
   * @return The best Solution attending to the comparator or <code>null<code>
   * if the SolutionSet is empty
   */
  public Solution best(Comparator comparator){
    int indexBest = indexBest(comparator);
    if (indexBest < 0) {
        return null;
    } else {
        return solutionsList_.get(indexBest);
    }
        
  } // best  

  
  /** 
   * Returns the index of the worst Solution using a <code>Comparator</code>.
   * If there are more than one occurrences, only the index of the first one is returned
   * @param comparator <code>Comparator</code> used to compare solutions.
   * @return The index of the worst Solution attending to the comparator or 
   * <code>-1<code> if the SolutionSet is empty
   */
  public int indexWorst(Comparator comparator){
    
   
    if ((solutionsList_ == null) || (this.solutionsList_.isEmpty())) {
        return -1;
    }
    
    int index = 0;
    Solution worstKnown = solutionsList_.get(0), candidateSolution;
    int flag;
    for (int i = 1; i < solutionsList_.size(); i++) {        
        candidateSolution = solutionsList_.get(i);
        flag = comparator.compare(worstKnown, candidateSolution);
        if (flag == -1) {
            index = i;
            worstKnown = candidateSolution;
        }
    }
    
    return index;
        
  } // indexWorst
  
  /** 
   * Returns the worst Solution using a <code>Comparator</code>.
   * If there are more than one occurrences, only the first one is returned
   * @param comparator <code>Comparator</code> used to compare solutions.
   * @return The worst Solution attending to the comparator or <code>null<code>
   * if the SolutionSet is empty
   */
  public Solution worst(Comparator comparator){
    
    int index = indexWorst(comparator);
    if (index < 0) {
        return null;
    } else {
        return solutionsList_.get(index);
    }
        
  } // worst
  
  
  /** 
   * Returns the number of solutions in the SolutionSet.
   * @return The size of the SolutionSet.
   */  
  public int size(){
    return solutionsList_.size();
  } // size
  
  public void printToFile(String path, String contents) throws JMException{
      
      try{
          FileOutputStream fos   = new FileOutputStream(path,true);
          OutputStreamWriter osw = new OutputStreamWriter(fos)    ;
          BufferedWriter bw      = new BufferedWriter(osw)        ;

          bw.write(contents);
          bw.newLine();
          
          /* Close the file */
          bw.close();
        }
        catch(IOException e){
          Configuration.logger_.severe("Error acceding to the file");
          e.printStackTrace();      
        }        
      
  }
  
  /**
   * Writes generations into set in a file
   * @param path The output file name
   * @param frontFormat The format of the output
   *        0 ... output only the objective vectors
   *        1 ... output decision and objective vectors
   */
  public void printGenerationsToFile(String path, int frontFormat) throws JMException{
   String s;
   
   try{
      FileOutputStream fos   = new FileOutputStream(path,true);
      OutputStreamWriter osw = new OutputStreamWriter(fos)    ;
      BufferedWriter bw      = new BufferedWriter(osw)        ;
          
      for (int i = 0; i < solutionsList_.size(); i++) {
          bw.write(i+1+ "\t");  //Sequence number
          if(frontFormat==1){
              int numberOfVariables = solutionsList_.get(0).getDecisionVariables().length;
              for (int j = 0; j < numberOfVariables; j++){  //Decision variables values [x1 x2 x3 ... xN]
                s=String.format("% f\t",solutionsList_.get(i).getDecisionVariables()[j].getValue());  
                bw.write(s);
              }   
          }
          bw.write(solutionsList_.get(i).toString()); //Objectives [f1 f2 ... fN]
          bw.newLine();
      }
      
      /* Close the file */
      bw.close();
    }
    catch(IOException e){
      Configuration.logger_.severe("Error acceding to the file");
      e.printStackTrace();      
    }    
  } //printGenerationsToFile
  
  
  /**
   * Writes the parameters of Problem, Algorithm and Output into set in a file
   * @param path The output file name
   * @param properties The parameters of the Problem, Algorithm and Output
   */
  public void printParametersToFile(String path, Properties properties, boolean mode){
    String key, value;  
      
    try{
      FileOutputStream fos   = new FileOutputStream(path,mode);
      OutputStreamWriter osw = new OutputStreamWriter(fos)    ;
      BufferedWriter bw      = new BufferedWriter(osw)        ;
      
      Set parametersSet = properties.keySet();
      Object [] keySet = parametersSet.toArray();
      Arrays.sort(keySet);
      
      for(int i=0; i<keySet.length; i++){
          key=(String) keySet[i];
          value = properties.getProperty(key);   
          bw.write("# " + key + " " + value);
          bw.newLine();
      }
        
      bw.newLine();
      /* Close the file */
      bw.close();
    }
    catch(IOException e){
      Configuration.logger_.severe("Error acceding to the file");
      e.printStackTrace();      
    }
  } //printParametersToFile
  
  /**
   * Writes the sequence number, decision variable values and objective function values
   * of the <code>Solution</code> objects into the set in a file
   * @param path The output file name
   */
  
  public void printVariablesAndObjectivesToFile(String path) throws JMException{
    String s;  
    try{
      FileOutputStream fos   = new FileOutputStream(path,true);
      OutputStreamWriter osw = new OutputStreamWriter(fos)    ;
      BufferedWriter bw      = new BufferedWriter(osw)        ;
      
      if(solutionsList_.size() == 0)
      	return;
      
      int numberOfVariables = solutionsList_.get(0).getDecisionVariables().length;
      int numberOfFeatures = Params.getNumberOfFeatures();
      
      for (int i = 0; i < solutionsList_.size(); i++) {  
          bw.write(i+1+ "\t");  //Sequence number
          
          for (int j = 0; j < numberOfVariables; j++){  //Decision variables values [x1 x2 x3 ... xN]
              s=String.format("% f\t",solutionsList_.get(i).getDecisionVariables()[j].getValue());
              bw.write(s);
          }
          
          for (int j = 0; j < Params.getNumberOfObjectives(); j++){
          //s= String.format("% d\t",solutionsList_.get(i).getObjective(j));
          s= String.valueOf(solutionsList_.get(i).getObjective(j)) + "\t";
          bw.write(s); //Objectives [f1 f2 ... fN]
          }
          
          for (int j = 0; j < numberOfFeatures; j++){  //Features [x1 x2 x3 ... xN]
              s=String.format("% f\t",solutionsList_.get(i).getFeature(j));
              bw.write(s);
          }
          
          bw.newLine();
      }
      bw.newLine();
      /* Close the file */
      bw.close();
    }
    catch(IOException e){
      Configuration.logger_.severe("Error acceding to the file");
      e.printStackTrace();      
    }
  } //printVariablesAndObjectsToFile
  
  /**
   * Writes the sequence number, decision variable values and objective function values
   * of the <code>Solution</code> objects into the set in a file
   * @param path The output file name
   * @param numberOfEvaluationsSoFar Sequence number of the first evaluation that will be printed
   */
  public void printVariablesAndObjectivesToFileFromSpecificSequenceNumber(String path, int numberOfEvaluationsSoFar) throws JMException{
	    String s;  
	    try{
	      FileOutputStream fos   = new FileOutputStream(path,true);
	      OutputStreamWriter osw = new OutputStreamWriter(fos)    ;
	      BufferedWriter bw      = new BufferedWriter(osw)        ;
	      
	      int numberOfVariables = solutionsList_.get(0).getDecisionVariables().length;
	      for (int i = 0; i < solutionsList_.size(); i++) {  
	          bw.write((i + numberOfEvaluationsSoFar) + "\t");  //Sequence number
	          
	          for (int j = 0; j < numberOfVariables; j++){  //Decision variables values [x1 x2 x3 ... xN]
	              s=String.format("% f\t",solutionsList_.get(i).getDecisionVariables()[j].getValue());
	              bw.write(s);
	          }
	          bw.write(solutionsList_.get(i).toString()); //Objectives [f1 f2 ... fN]
	          bw.newLine();
	      }
	      /* Close the file */
	      bw.close();
	    }
	    catch(IOException e){
	      Configuration.logger_.severe("Error acceding to the file");
	      e.printStackTrace();      
	    }
	  } //printVariablesAndObjectsToFile
  
  /** 
   * Writes the objective function values of the <code>Solution</code> 
   * objects into the set in a file.
   * @param path The output file name
   */ 
  public void printObjectivesToFile(String path){
    try {
      /* Open the file */
      FileOutputStream fos   = new FileOutputStream(path)     ;
      OutputStreamWriter osw = new OutputStreamWriter(fos)    ;
      BufferedWriter bw      = new BufferedWriter(osw)        ;
                        
      for (int i = 0; i < solutionsList_.size(); i++) {
        //if (this.vector[i].getFitness()<1.0) {
        bw.write(solutionsList_.get(i).toString());
        bw.newLine();
        //}
      }
      
      /* Close the file */
      bw.close();
    }catch (IOException e) {
      Configuration.logger_.severe("Error acceding to the file");
      e.printStackTrace();
    }
  } // printObjectivesToFile
 
  /**
   * Writes the decision variable values of the <code>Solution</code>
   * solutions objects into the set in a file.
   * @param path The output file name
   */
  public void printVariablesToFile(String path){
    try {
      /* Open the file */
      FileOutputStream fos   = new FileOutputStream(path)     ;
      OutputStreamWriter osw = new OutputStreamWriter(fos)    ;
      BufferedWriter bw      = new BufferedWriter(osw)        ;            
            
      int numberOfVariables = solutionsList_.get(0).getDecisionVariables().length ;
      for (int i = 0; i < solutionsList_.size(); i++) {  
      	for (int j = 0; j < numberOfVariables; j++)
          bw.write(solutionsList_.get(i).getDecisionVariables()[j].toString() + " ");
        bw.newLine();        
      }
      
      /* Close the file */
      bw.close();
    }catch (IOException e) {
      Configuration.logger_.severe("Error acceding to the file");
      e.printStackTrace();
    }       
  } // printVariablesToFile
     
  /** 
   * Empties the SolutionSet
   */
  public void clear(){
    solutionsList_.clear();
  } // clear
    
  /** 
   * Deletes the <code>Solution</code> at position i in the set.
   * @param i The position of the solution to remove.
   */
  public void remove(int i){        
    if (i > solutionsList_.size()-1) {            
      Configuration.logger_.severe("Size is: "+this.size());
    } // if
    solutionsList_.remove(i);    
  } // remove
    
    
  /**
   * Returns an <code>Iterator</code> to access to the solution set list.
   * @return the <code>Iterator</code>.
   */    
  public Iterator<Solution> iterator(){
    return solutionsList_.iterator();
  } // iterator   
   
  /** 
   * Returns a new <code>SolutionSet</code> which is the result of the union
   * between the current solution set and the one passed as a parameter.
   * @param solutionSet SolutionSet to join with the current solutionSet.
   * @return The result of the union operation.
   */
  public SolutionSet union(SolutionSet solutionSet) {       
    //Check the correct size. In development
    int newSize = this.size() + solutionSet.size();
    if (newSize < capacity_)
      newSize = capacity_;
        
    //Create a new population 
    SolutionSet union = new SolutionSet(newSize);                
    for (int i = 0; i < this.size(); i++) {      
      union.add(this.get(i));
    } // for
        
    for (int i = this.size(); i < (this.size() + solutionSet.size()); i++) {
      union.add(solutionSet.get(i-this.size()));
    } // for
            
    return union;        
  } // union                   
    
  /** 
   * Replaces a solution by a new one
   * @param position The position of the solution to replace
   * @param solution The new solution
   */
  public void replace(int position, Solution solution) {
    if (position > this.solutionsList_.size()) {
      solutionsList_.add(solution);
    } // if 
    solutionsList_.remove(position);
    solutionsList_.add(position,solution);
  } // replace

  /**
   * Copies the objectives of the solution set to a matrix
   * @return A matrix containing the objectives
   */
  public double [][] writeObjectivesToMatrix() {
    if (this.size() == 0) {
      return null;
    }
    double [][] objectives;
    objectives = new double[size()][get(0).numberOfObjectives()];
    for (int i = 0; i < size(); i++) {
      for (int j = 0; j < get(0).numberOfObjectives(); j++) {
        objectives[i][j] = get(i).getObjective(j);
      }
    }
    return objectives;
  } // writeObjectivesMatrix

  /**
   * Prints the info about generation in a file
   * @param path The output file name
   * @param problem The Problem instance
   * @param population Current population of solutions
   * @param populationSize size of the population
   * @param iteration Iteration number of generation
   * @param numberOfChildBetter How many times was the child better than the parent
   * @param numberOfParentBetter How many times was the parent better than the child
   * @param numberOfIncomparable How many times was the child and the parent incomparable
   * @param hiperSwitch If 1 calculate hypervolume, if 0 dont calculate it
   */
public void printInfoAboutGenerations(String path, Problem problem, SolutionSet population, int populationSize, int iteration, int numberOfChildBetter, int numberOfParentBetter, int numberOfIncomparable, int additionalEvaluations, double NNAccuracy) {
	String contents;
	int numberOfFeasible = 0;
	int numberOfInfeasible = 0;
	int numberOfNonCalculated =0;
	double solutionViolation;
	double constraintViolation = 0;
	double hypervolume;
	double spread;
	
	try{
        FileOutputStream fos   = new FileOutputStream(path,true);
        OutputStreamWriter osw = new OutputStreamWriter(fos)    ;
        BufferedWriter bw      = new BufferedWriter(osw)        ;

        contents=String.format("%d\t", iteration);
        
        for(int i=0; i< populationSize; i++)
        {
        	solutionViolation = population.get(i).getOverallConstraintViolation();
        	
        	if(Math.abs(solutionViolation) > 0) 
            {
            	numberOfInfeasible++;
            	constraintViolation = constraintViolation + Math.abs(solutionViolation);
            }
            else
            	numberOfFeasible++;
        }
        
        contents=contents + String.format("%d\t", numberOfFeasible);
        contents=contents + String.format("%d\t", numberOfInfeasible);
        contents=contents + String.format("%d\t", numberOfNonCalculated);
        if(numberOfInfeasible == 0)
      	  contents=contents + String.format("%d\t", 0);
        else
      	  contents=contents + String.format("%s\t", constraintViolation/numberOfInfeasible);
        
        Ranking ranking = new Ranking(population);
        contents=contents + String.format("%d\t", ranking.getSubfront(0).size());
        double [][] front = ranking.getSubfront(0).writeObjectivesToMatrix();
        
        if(Params.getHypervolumeSwitch() == 1) { 
      	  Hypervolume hyp = new Hypervolume();
    
      	  //hypervolume = hyp.calculateHypervolume(front, ranking.getSubfront(0).size(), problem.numberOfObjectives_);
      	  hypervolume = hyp.hypervolume(front, Params.getFrontForHypervolume(), problem.numberOfObjectives_);

        }
        else
      	  hypervolume = 0;
        
        contents = contents + String.format("%s\t", hypervolume);
        
        GeneralizedSpread spr = new GeneralizedSpread();
        spread = spr.generalizedSpread(front, Params.getFrontForHypervolume(), problem.numberOfObjectives_);
        contents = contents + String.format("%s\t", spread);        
        
        contents=contents + String.format("%d\t", numberOfChildBetter);
        contents=contents + String.format("%d\t", numberOfParentBetter);
        contents=contents + String.format("%d\t", numberOfIncomparable);
        contents=contents + String.format("%d\t", additionalEvaluations);
        if(NNAccuracy != -1)
        	contents=contents + String.format("%s\t", NNAccuracy);
        else
        	contents=contents + String.format("%s\t", "");
        
        bw.write(contents);
        bw.newLine();
        
        /* Close the file */
        bw.close();
      }
      catch(IOException e){
        Configuration.logger_.severe("Error acceding to the file");
        e.printStackTrace();      
      }          
	
} //printInfoAboutGenerations

public void printInfoAboutGenerations(String path, Problem problem, SolutionSet population, int populationSize, int iteration, int numberOfChildBetter, int numberOfParentBetter, int numberOfIncomparable, int additionalEvaluations, int exactEvaluationsDuringEvolutionProcess) {
	String contents;
	int numberOfFeasible = 0;
	int numberOfInfeasible = 0;
	int numberOfNonCalculated =0;
	double solutionViolation;
	double constraintViolation = 0;
	double hypervolume;
	double spread;
	
	try{
        FileOutputStream fos   = new FileOutputStream(path,true);
        OutputStreamWriter osw = new OutputStreamWriter(fos)    ;
        BufferedWriter bw      = new BufferedWriter(osw)        ;

        contents=String.format("%d\t", iteration);
        
        for(int i=0; i< populationSize; i++)
        {
        	solutionViolation = population.get(i).getOverallConstraintViolation();
        	
        	if(Math.abs(solutionViolation) > 0) 
            {
            	numberOfInfeasible++;
            	constraintViolation = constraintViolation + Math.abs(solutionViolation);
            }
            else
            	numberOfFeasible++;
        }
        
        contents=contents + String.format("%d\t", numberOfFeasible);
        contents=contents + String.format("%d\t", numberOfInfeasible);
        contents=contents + String.format("%d\t", numberOfNonCalculated);
        if(numberOfInfeasible == 0)
      	  contents=contents + String.format("%d\t", 0);
        else
      	  contents=contents + String.format("%s\t", constraintViolation/numberOfInfeasible);
        
        Ranking ranking = new Ranking(population);
        contents=contents + String.format("%d\t", ranking.getSubfront(0).size());
        double [][] front = ranking.getSubfront(0).writeObjectivesToMatrix();
        
        if(Params.getHypervolumeSwitch() == 1) { 
      	  Hypervolume hyp = new Hypervolume();

      	  //hypervolume = hyp.calculateHypervolume(front, ranking.getSubfront(0).size(), problem.numberOfObjectives_);  m    	          
      	  
      	  hypervolume = hyp.hypervolume(front, Params.getFrontForHypervolume(), problem.numberOfObjectives_);
        }
        else
      	  hypervolume = 0;
        
        contents = contents + String.format("%s\t", hypervolume);
        
        GeneralizedSpread spr = new GeneralizedSpread();
        spread = spr.generalizedSpread(front, Params.getFrontForHypervolume(), problem.numberOfObjectives_);
        contents = contents + String.format("%s\t", spread); 
        
        contents=contents + String.format("%d\t", numberOfChildBetter);
        contents=contents + String.format("%d\t", numberOfParentBetter);
        contents=contents + String.format("%d\t", numberOfIncomparable);
        contents=contents + String.format("%d\t", additionalEvaluations);
        contents=contents + String.format("%d\t", exactEvaluationsDuringEvolutionProcess);     
        contents=contents + String.format("%s\t", getCurrentExecutionTime(System.currentTimeMillis()));
        
        bw.write(contents);
        bw.newLine();
        
        /* Close the file */
        bw.close();
      }
      catch(IOException e){
        Configuration.logger_.severe("Error acceding to the file");
        e.printStackTrace();      
      }          
	
} //printInfoAboutGenerations

private String getCurrentExecutionTime(long currentTimeMillis) {
	
	long initTime = Params.getStartTime();
  long diff = currentTimeMillis - initTime;
  
  long secondInMillis = 1000;
  long minuteInMillis = secondInMillis * 60;
  long hourInMillis = minuteInMillis * 60;
  long dayInMillis = hourInMillis * 24;
  long yearInMillis = dayInMillis * 365;

  String elapsedYears = "0" + String.valueOf(diff / yearInMillis);
  diff = diff % yearInMillis;
  String elapsedDays = "0" + String.valueOf(diff / dayInMillis);
  diff = diff % dayInMillis;
  String elapsedHours = "0" + String.valueOf(diff / hourInMillis);
  diff = diff % hourInMillis;
  String elapsedMinutes = "0" + String.valueOf(diff / minuteInMillis);
  diff = diff % minuteInMillis;
  String elapsedSeconds = "0" + String.valueOf(diff / secondInMillis);

  String evaluationTime = elapsedHours.substring(elapsedHours.length()-2)  +":" + elapsedMinutes.substring(elapsedMinutes.length()-2) + ":" + elapsedSeconds.substring(elapsedSeconds.length()-2);
  
	
	return evaluationTime;
}
} // SolutionSet

