//  Solution.java
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>
//
//  Description: 
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

import jmetal.encodings.variable.Binary;
import jmetal.util.Configuration;
import jmetal.util.Configuration.*;

/**
 * Class representing a solution for a problem.
 */
public class Solution implements Serializable {  
	/**
	 * Stores the problem 
	 */
	Problem problem_ ;
	
  /**
   * Stores the type of the variable
   */	
  private SolutionType type_ ; 

  /**
   * Stores the decision variables of the solution.
   */
  private Variable[] variable_ ;

  /**
   * Stores the objectives values of the solution.
   */
  private double [] objective_ ;
  
  /**
   * Stores the values of features of the solution.
   */
  private double [] features_ ;
  
  /**
   * Stores the values of variance of objectives of the solution.
   */
  private double [] standardDeviance_ ;

  /**
   * Stores the number of objective values of the solution
   */
  private int numberOfObjectives_ ;

  /**
   * Stores the number of feature values of the solution
   */
  private int numberOfFeatures_ ;
  
  /**
   * Stores the so called fitness value. Used in some metaheuristics
   */
  private double fitness_ ;

  /**
   * Used in algorithm AbYSS, this field is intended to be used to know
   * when a <code>Solution</code> is marked.
   */
  private boolean marked_ ;
  
  /**
   * Used with meta-models when a <code>Solution</code> is approximated.
   */
  private boolean exactlyEvaluated_ ;
  

  /**
   * Stores the so called rank of the solution. Used in NSGA-II
   */
  private int rank_ ;

  /**
   * Stores the overall constraint violation of the solution.
   */
  private double  overallConstraintViolation_ ;

  /**
   * Stores the number of constraints violated by the solution.
   */
  private int  numberOfViolatedConstraints_ ;

  /**
   * Stores the solution feasibility. Value -1 = (probably) infeasible. 0 = uncertainly feasible.   1 = (probably) feasible 
   */
  private int  solutionFeasibility_ ;
  
  /**
   * This field is intended to be used to know the location of
   * a solution into a <code>SolutionSet</code>. Used in MOCell
   */
  private int location_ ;

  /**
   * Stores the distance to his k-nearest neighbor into a 
   * <code>SolutionSet</code>. Used in SPEA2.
   */
  private double kDistance_ ; 

  /**
   * Stores the crowding distance of the the solution in a 
   * <code>SolutionSet</code>. Used in NSGA-II.
   */
  private double crowdingDistance_ ; 

  /**
   * Stores the distance between this solution and a <code>SolutionSet</code>.
   * Used in AbySS.
   */
  private double distanceToSolutionSet_ ;       

  /**
   * Constructor.
   */
  public Solution() {        
    problem_                      = null  ;
    marked_                       = false ;
    exactlyEvaluated_			  = true;
    overallConstraintViolation_   = 0.0   ;
    numberOfViolatedConstraints_  = 0     ;  
    type_                         = null ;
    variable_                     = null ;
    objective_                    = null ;
    features_					  					= null ;
    standardDeviance_							= null ;
  } // Solution

  /**
   * Constructor
   * @param numberOfObjectives Number of objectives of the solution
   * 
   * This constructor is used mainly to read objective values from a file to
   * variables of a SolutionSet to apply quality indicators
   */
  public Solution(int numberOfObjectives) {
    numberOfObjectives_ = numberOfObjectives;
    objective_          = new double[numberOfObjectives];
    features_           = new double[numberOfObjectives];
    standardDeviance_   = new double[numberOfObjectives];
  }
  
  /** 
   * Constructor.
   * @param problem The problem to solve
   * @throws ClassNotFoundException 
   */
  public Solution(Problem problem) throws ClassNotFoundException{
    problem_ = problem ; 
    type_ = problem.getSolutionType() ;
    numberOfObjectives_ = problem.getNumberOfObjectives() ;
    numberOfFeatures_ = problem.getNumberOfFeatures();
    objective_          = new double[numberOfObjectives_] ;
    features_          = new double[numberOfFeatures_] ;
    standardDeviance_  = new double[numberOfObjectives_] ;
    solutionFeasibility_ =1; 		//by Miha. Because if not calculated (conf. int.) differently, we assume that it is ok. There is still no of violated constraints for infeasibility

    // Setting initial values
    fitness_              = 0.0 ;
    kDistance_            = 0.0 ;
    crowdingDistance_     = 0.0 ;        
    distanceToSolutionSet_ = Double.POSITIVE_INFINITY ;
    //<-

    //variable_ = problem.solutionType_.createVariables() ; 
    variable_ = type_.createVariables() ; 
  } // Solution
  
  static public Solution getNewSolution(Problem problem) throws ClassNotFoundException {
    return new Solution(problem) ;
  }
  
  /** 
   * Constructor
   * @param problem The problem to solve
   */
  public Solution(Problem problem, Variable [] variables){
    problem_ = problem ;
    type_ = problem.getSolutionType() ;
    numberOfObjectives_ = problem.getNumberOfObjectives() ;
    numberOfFeatures_ = problem.getNumberOfFeatures();
    objective_          = new double[numberOfObjectives_] ;
    features_          = new double[numberOfFeatures_] ;
    standardDeviance_   = new double[numberOfObjectives_] ;

    // Setting initial values
    fitness_              = 0.0 ;
    kDistance_            = 0.0 ;
    crowdingDistance_     = 0.0 ;        
    distanceToSolutionSet_ = Double.POSITIVE_INFINITY ;
    //<-

    //variable_ = variables ; jMetal
    variable_ = variables.clone(); //Rok Prodan modifications 
  } // Constructor
  
  /** 
   * Copy constructor.
   * @param solution Solution to copy.
   */    
  public Solution(Solution solution) {            
    problem_ = solution.problem_ ;
    type_ = solution.type_;

    numberOfObjectives_ = solution.numberOfObjectives();    
    objective_ = new double[numberOfObjectives_];
    for (int i = 0; i < objective_.length;i++) {
      objective_[i] = solution.getObjective(i);
    } // for
    //<-
    numberOfFeatures_ = solution.numberOfFeatures();
    features_ = new double[numberOfFeatures_];
    for (int i = 0; i < features_.length;i++) {
    	features_[i] = solution.getFeature(i);
    } // for
    //<-       
    standardDeviance_ = new double[numberOfObjectives_];
    for (int i = 0; i < standardDeviance_.length;i++) {
    	standardDeviance_[i] = solution.getstandardDeviance(i);
    } // for
    //<-
    
    
    
    variable_ = type_.copyVariables(solution.variable_) ;
    overallConstraintViolation_  = solution.getOverallConstraintViolation();
    numberOfViolatedConstraints_ = solution.getNumberOfViolatedConstraint();
    distanceToSolutionSet_ = solution.getDistanceToSolutionSet();
    crowdingDistance_     = solution.getCrowdingDistance();
    kDistance_            = solution.getKDistance();                
    fitness_              = solution.getFitness();
    marked_               = solution.isMarked();
    exactlyEvaluated_     = solution.isExactllyEvaluated();
    rank_                 = solution.getRank();
    location_             = solution.getLocation();
    solutionFeasibility_  = solution.getSolutionFeasibility();
  } // Solution

  /**
   * Sets the distance between this solution and a <code>SolutionSet</code>.
   * The value is stored in <code>distanceToSolutionSet_</code>.
   * @param distance The distance to a solutionSet.
   */
  public void setDistanceToSolutionSet(double distance){
    distanceToSolutionSet_ = distance;
  } // SetDistanceToSolutionSet

  /**
   * Gets the distance from the solution to a <code>SolutionSet</code>. 
   * <b> REQUIRE </b>: this method has to be invoked after calling 
   * <code>setDistanceToPopulation</code>.
   * @return the distance to a specific solutionSet.
   */
  public double getDistanceToSolutionSet(){
    return distanceToSolutionSet_;
  } // getDistanceToSolutionSet


  /** 
   * Sets the distance between the solution and its k-nearest neighbor in 
   * a <code>SolutionSet</code>. The value is stored in <code>kDistance_</code>.
   * @param distance The distance to the k-nearest neighbor.
   */
  public void setKDistance(double distance){
    kDistance_ = distance;
  } // setKDistance

  /** 
   * Gets the distance from the solution to his k-nearest nighbor in a 
   * <code>SolutionSet</code>. Returns the value stored in
   * <code>kDistance_</code>. <b> REQUIRE </b>: this method has to be invoked 
   * after calling <code>setKDistance</code>.
   * @return the distance to k-nearest neighbor.
   */
  public double getKDistance(){
    return kDistance_;
  } // getKDistance

  /**
   * Sets the crowding distance of a solution in a <code>SolutionSet</code>.
   * The value is stored in <code>crowdingDistance_</code>.
   * @param distance The crowding distance of the solution.
   */  
  public void setCrowdingDistance(double distance){
    crowdingDistance_ = distance;
  } // setCrowdingDistance


  /** 
   * Gets the crowding distance of the solution into a <code>SolutionSet</code>.
   * Returns the value stored in <code>crowdingDistance_</code>.
   * <b> REQUIRE </b>: this method has to be invoked after calling 
   * <code>setCrowdingDistance</code>.
   * @return the distance crowding distance of the solution.
   */
  public double getCrowdingDistance(){
    return crowdingDistance_;
  } // getCrowdingDistance

  /**
   * Sets the fitness of a solution.
   * The value is stored in <code>fitness_</code>.
   * @param fitness The fitness of the solution.
   */
  public void setFitness(double fitness) {
    fitness_ = fitness;
  } // setFitness

  /**
   * Gets the fitness of the solution.
   * Returns the value of stored in the variable <code>fitness_</code>.
   * <b> REQUIRE </b>: This method has to be invoked after calling 
   * <code>setFitness()</code>.
   * @return the fitness.
   */
  public double getFitness() {
    return fitness_;
  } // getFitness

  /**
   * Sets the value of the i-th objective.
   * @param i The number identifying the objective.
   * @param value The value to be stored.
   */
  public void setObjective(int i, double value) {
    objective_[i] = value;
  } // setObjective

  /**
   * Returns the value of the i-th objective.
   * @param i The value of the objective.
   */
  public double getObjective(int i) {
    return objective_[i];
  } // getObjective

  /**
   * Returns the number of objectives.
   * @return The number of objectives.
   */
  public int numberOfObjectives() {
    if (objective_ == null)
      return 0 ;
    else
      return numberOfObjectives_;
  } // numberOfObjectives

  /**
   * Sets the value of the i-th variance.
   * @param i The number identifying the objective.
   * @param value The variance value to be stored.
   */
  public void setStandardDeviance(int i, double value) {
  	standardDeviance_[i] = value;
  } // setVariance

  /**
   * Sets the value of the i-th variance.
   * @param i The number identifying the objective.
   * @param value The variance value to be stored.
   */
  public void setStandardDevianceToZero() {
    for(int i=0; i<numberOfObjectives_; i++)
    	standardDeviance_[i] = 0;
  } // setVarianceToZero
  
  /**
   * Returns the value of the i-th objective variance.
   * @param i The value of the objective variance.
   */
  public double getstandardDeviance(int i) {
    return standardDeviance_[i];
  } // getstandardDeviance

  /**
   * Returns the value of the i-th objective variance.
   * @param i The value of the objective variance.
   */
  public double[] getstandardDeviances() {
    return standardDeviance_;
  } // getstandardDeviances
  
  /**
   * Sets the value of the i-th feature.
   * @param i The number identifying the feature.
   * @param value The value to be stored.
   */
  public void setFeature(int i, double value) {
    features_[i] = value;
  } // setFeature
  
  /**
   * Returns the value of the i-th feature.
   * @param i The value of the feature.
   */
  public double getFeature(int i) {
    return features_[i];
  } // getFeature
  
  /**
   * Returns the number of features.
   * @return The number of features.
   */
  public int numberOfFeatures() {
    if (features_ == null)
      return 0 ;
    else
      return numberOfFeatures_;
  } // numberOfObjectives
  
  /**  
   * Returns the number of decision variables of the solution.
   * @return The number of decision variables.
   */
  public int numberOfVariables() {
    return problem_.getNumberOfVariables() ;
  } // numberOfVariables

  /** 
   * Returns a string representing the solution.
   * @return The string.
   */
  public String toString() {
    String s;  
    String aux="";
    for (int i = 0; i < this.numberOfObjectives_; i++){
      s=String.format("% f\t", this.getObjective(i));
      aux = aux + s;
    }  

    return aux;
  } // toString

  /**
   * Returns the decision variables of the solution.
   * @return the <code>DecisionVariables</code> object representing the decision
   * variables of the solution.
   */
  public Variable[] getDecisionVariables() {
    return variable_ ;
  } // getDecisionVariables

  /**
   * Sets the decision variables for the solution.
   * @param decisionVariables The <code>DecisionVariables</code> object 
   * representing the decision variables of the solution.
   */
  public void setDecisionVariables(Variable [] variables) {
    variable_ = variables ;
  } // setDecisionVariables

  /**
   * Indicates if the solution is marked.
   * @return true if the method <code>marked</code> has been called and, after 
   * that, the method <code>unmarked</code> hasn't been called. False in other
   * case.
   */
  public boolean isMarked() {
    return this.marked_;
  } // isMarked

  /**
   * Establishes the solution as marked.
   */
  public void marked() {
    this.marked_ = true;
  } // marked

  /**
   * Established the solution as unmarked.
   */
  public void unMarked() {
    this.marked_ = false;
  } // unMarked

  
  public boolean isExactllyEvaluated() {
	    return this.exactlyEvaluated_;
	  } // isExactllyEvaluated
  
  
  public void setExactllyEvaluated(boolean b) {
    this.exactlyEvaluated_ = b;
    
  } // marked

  
  /**  
   * Sets the rank of a solution. 
   * @param value The rank of the solution.
   */
  public void setRank(int value){
    this.rank_ = value;
  } // setRank

  /**
   * Gets the rank of the solution.
   * <b> REQUIRE </b>: This method has to be invoked after calling 
   * <code>setRank()</code>.
   * @return the rank of the solution.
   */
  public int getRank(){
    return this.rank_;
  } // getRank

  /**
   * Sets the overall constraints violated by the solution.
   * @param value The overall constraints violated by the solution.
   */
  public void setOverallConstraintViolation(double value) {
    this.overallConstraintViolation_ = value;
  } // setOverallConstraintViolation

  /**
   * Gets the overall constraint violated by the solution.
   * <b> REQUIRE </b>: This method has to be invoked after calling 
   * <code>overallConstraintViolation</code>.
   * @return the overall constraint violation by the solution.
   */
  public double getOverallConstraintViolation() {
    return this.overallConstraintViolation_;
  }  //getOverallConstraintViolation


  /**
   * Sets the number of constraints violated by the solution.
   * @param value The number of constraints violated by the solution.
   */
  public void setNumberOfViolatedConstraint(int value) {
    this.numberOfViolatedConstraints_ = value;
  } //setNumberOfViolatedConstraint

  /**
   * Gets the number of constraint violated by the solution.
   * <b> REQUIRE </b>: This method has to be invoked after calling
   * <code>setNumberOfViolatedConstraint</code>.
   * @return the number of constraints violated by the solution.
   */
  public int getNumberOfViolatedConstraint() {
    return this.numberOfViolatedConstraints_;
  } // getNumberOfViolatedConstraint

  
  /**
   * Sets the knownFeasibility of the solution.
   * @param value -1 = (probably) infeasible. 0 = uncertainly feasible.   1 = (probably) feasible 
   */
  public void setSolutionFeasibility(int value) {
    this.solutionFeasibility_ = value;
  } //setUnknownFeasibility

  /**
   * Gets the knownFeasibility of the solution.
   * <b> REQUIRE </b>: This method has to be invoked after calling
   * <code>setUnknownFeasibility</code>.
   * @return value -1 = (probably) infeasible. 0 = uncertainly feasible.   1 = (probably) feasible 
   */
  public int getSolutionFeasibility() {
    return this.solutionFeasibility_;
  } // getUnknownFeasibility
  
  
  /**
   * Sets the location of the solution into a solutionSet. 
   * @param location The location of the solution.
   */
  public void setLocation(int location) {
    this.location_ = location;
  } // setLocation

  /**
   * Gets the location of this solution in a <code>SolutionSet</code>.
   * <b> REQUIRE </b>: This method has to be invoked after calling
   * <code>setLocation</code>.
   * @return the location of the solution into a solutionSet
   */
  public int getLocation() {
    return this.location_;
  } // getLocation

  /**
   * Sets the type of the variable. 
   * @param type The type of the variable.
   */
  //public void setType(String type) {
   // type_ = Class.forName("") ;
  //} // setType

  /**
   * Sets the type of the variable. 
   * @param type The type of the variable.
   */
  public void setType(SolutionType type) {
    type_ = type ;
  } // setType

  /**
   * Gets the type of the variable
   * @return the type of the variable
   */
  public SolutionType getType() {
    return type_;
  } // getType

  /** 
   * Returns the aggregative value of the solution
   * @return The aggregative value.
   */
  public double getAggregativeValue() {
    double value = 0.0;                
    for (int i = 0; i < numberOfObjectives(); i++){            
      value += getObjective(i);
    }                
    return value;
  } // getAggregativeValue

  /**
   * Returns the number of bits of the chromosome in case of using a binary
   * representation
   * @return The number of bits if the case of binary variables, 0 otherwise
   */
  public int getNumberOfBits() {
    int bits = 0 ;
    
    for (int i = 0;  i < variable_.length  ; i++)
	    try {
	      if ((variable_[i].getVariableType() == Class.forName("jmetal.base.variable.Binary")) ||
	          (variable_[i].getVariableType() == Class.forName("jmetal.base.variable.BinaryReal")))
	        bits += ((Binary)(variable_[i])).getNumberOfBits() ;
      } catch (ClassNotFoundException e) {
	      e.printStackTrace();
      }
    
    return bits ;
  } // getNumberOfBits
  
  /**
   * Writes decision variables values and objective values to file
   * @param path The output file name
   */
  public void printVariablesandObjectivesToFile(String path, int seqNumber){
  try{
      FileOutputStream fos   = new FileOutputStream(path,true);
      OutputStreamWriter osw = new OutputStreamWriter(fos)    ;
      BufferedWriter bw      = new BufferedWriter(osw)        ;
      
      bw.write(seqNumber + " | ");
      
      for (int i = 0; i < variable_.length; i++) { //Decision variables values [x1 x2 x3 ... xN] 
        bw.write(variable_[i] + " ");
      }
      
      bw.write(" | ");
      
      for (int i = 0; i < objective_.length; i++){ //Objectives [f1 f2 ... fN] 
        bw.write(objective_[i] + " ");
      }
      
      bw.write(" | ");
      
      for (int i = 0; i < features_.length; i++){ //Objectives [f1 f2 ... fN] 
        bw.write(features_[i] + " ");
      }
      
      bw.newLine();
      /* Close the file */
      bw.close();
    }
    catch(IOException e){
      Configuration.logger_.severe("Error acceding to the file");
      e.printStackTrace();      
    }    
           
  } //printVariablesandObjectivesToFile
  
  /**
   * Writes decision variables values and objective values and standard deviation to file
   * @param path The output file name
   */
  public void printVariablesandObjectivesAndStandardDeviationToFile(String path, int seqNumber, boolean exactlyEvaluated, double[] stdDeviances, String[] whichModelIsBetter, double [] solutionApproximatedValue, double[] solutionConfidenceInterval, String[] allObjectivesAndConfidences){
  try{
      FileOutputStream fos   = new FileOutputStream(path,true);
      OutputStreamWriter osw = new OutputStreamWriter(fos)    ;
      BufferedWriter bw      = new BufferedWriter(osw)        ;
      
      bw.write(seqNumber + " | ");
      
      for (int i = 0; i < variable_.length; i++) { //Decision variables values [x1 x2 x3 ... xN] 
        bw.write(variable_[i] + " ");
      }
      
      bw.write(" | ");
      
      for (int i = 0; i < objective_.length; i++){ //Objectives [f1 f2 ... fN] 
        bw.write(objective_[i] + " ");
      }
      
      if(features_.length >0)
      	bw.write(" | ");
      
      for (int i = 0; i < features_.length; i++){ //Objectives [f1 f2 ... fN] 
        bw.write(features_[i] + " ");
      }
            
      
      if(exactlyEvaluated == true)
      {
      	bw.write(" | Exact");      	
	      if(solutionApproximatedValue != null)
	      {
	      	bw.write(" | Aprox.obj.");
	      	bw.write(" | ");
	      	for(int i=0; i<solutionApproximatedValue.length; i++)
	      	{	      		
	      		bw.write(solutionApproximatedValue[i] + " ");
	      	}
	      	bw.write(" | Conf.inter.");
      		bw.write(" | ");
	      	for(int i=0; i<solutionConfidenceInterval.length; i++)
	      	{
	      		bw.write(solutionConfidenceInterval[i] + " ");
	      	}
	      	
	      	if(isSolutionInsideConfidenceInterval(solutionApproximatedValue, solutionConfidenceInterval))
	      		bw.write("--OK-- ");
	      	else
	      		bw.write("--WrongApproximation-- ");
	      
      	}
      }
      
      else
      {
      	bw.write(" | Aprox.");
      	//bw.write(" | ");
      	//for(int i=0; i<stdDeviances.length; i++)
      	//{      		
      	//	bw.write(stdDeviances[i] + " ");
      	//}
      	
      	if(solutionApproximatedValue != null)
	      {
	      	bw.write(" | Aprox.obj.");
	      	bw.write(" | ");
	      	for(int i=0; i<solutionApproximatedValue.length; i++)
	      	{	      		
	      		bw.write(solutionApproximatedValue[i] + " ");
	      	}
	      	bw.write(" | Conf.inter.");
      		bw.write(" | ");
	      	for(int i=0; i<solutionConfidenceInterval.length; i++)
	      	{
	      		bw.write(solutionConfidenceInterval[i] + " ");
	      	}
	      	
	      	bw.write("No_data ");
	      	
      	}
      	
      }
     
      
      if(whichModelIsBetter != null)
      	for(int i=0; i<objective_.length; i++)
      		bw.write(whichModelIsBetter[i] + " ");
      
      //Write detailed approximations
      if(whichModelIsBetter != null)
      	for(int i=0; i<allObjectivesAndConfidences.length;i++)
      		bw.write(allObjectivesAndConfidences[i]);
      
      bw.newLine();
      /* Close the file */
      bw.close();
    }
    catch(IOException e){
      Configuration.logger_.severe("Error acceding to the file");
      e.printStackTrace();      
    }    
           
  } //printVariablesandObjectivesAndStandardDeviationToFile

	private boolean isSolutionInsideConfidenceInterval(double[] solutionApproximatedValue, double[] solutionConfidenceInterval) {
	  
		boolean insideInterval = true;
		
		for(int i=0; i<objective_.length; i++)
		{
			if((objective_[i] < solutionApproximatedValue[i] - solutionConfidenceInterval[i]) || (objective_[i] > solutionApproximatedValue[i] + solutionConfidenceInterval[i]))
			{
				insideInterval = false;
				break;
			}
		}
		
	  return insideInterval;
  }
} // Solution
