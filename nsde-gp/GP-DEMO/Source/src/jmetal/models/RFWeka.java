//  GP.java
//
//  Author:
//       Dejan Petelin <dejan.petelin@ijs.si>
//
//  Copyright (c) 2011 Dejan Petelin
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

package jmetal.models;

import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.List;


import com.mathworks.toolbox.javabuilder.MWCellArray;
import com.mathworks.toolbox.javabuilder.MWException;
import com.mathworks.toolbox.javabuilder.MWNumericArray;

import sun.misc.BASE64Encoder;

import jmetal.core.*;
import jmetal.init.Params;
import jmetal.util.JMException;
import gpml.GPML;

import weka.classifiers.Classifier;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.classifiers.trees.RandomForest;
/**
 * Class representing model Random Forest
 */
public class RFWeka extends Model {  
  
	RandomForest[] tableOfRFClassifiers;
	int numberOfVariables_;
	int numberOfObjectives_;
	FastVector solutionParameters;
	Instances[] tableOfSolutionsForEachObjective;
	
	int capacity = 100;
	
  /** 
   * Constructor.
   * Creates a default instance of the GP model.
   * @throws Exception 
   */
  public RFWeka(int noInput, int noOutput) throws Exception {
  	super(noInput, noOutput);
  	  	
  	numberOfVariables_   = noInput;		
	  numberOfObjectives_  = noOutput;
	  
	  tableOfSolutionsForEachObjective = new Instances[numberOfObjectives_];
	  tableOfRFClassifiers = new RandomForest[numberOfObjectives_];
	  
	  initializeRandomForests();
	  InitializeTrainingSets();
  } // RF

	@Override
  public void update(SolutionSet solutions) {
		
		for(int i=0; i<solutions.size();i++)	    
      update(solutions.get(i));
	
  }

	@Override
  public void update(Solution solution) {
		
		try {
	    addSolutionToTrainingSet(solution);
    } 
		catch (JMException e) {
	    e.printStackTrace();
    }
		
		try {
	    buildClassifiers();
    }
		catch (Exception e) {	 
	    e.printStackTrace();
    }
  }

	private void buildClassifiers() throws Exception {

		for(int i=0; i< numberOfObjectives_; i++)
			tableOfRFClassifiers[i].buildClassifier(tableOfSolutionsForEachObjective[i]);		//Classifier cModel = (Classifier)new NaiveBayes();
			  
  }

	private void addSolutionToTrainingSet(Solution solution) throws JMException {
		
		 // Create the instance
		 Instance iExample;
		 
		 for(int i=0; i< numberOfObjectives_; i++)
		 {
			 iExample = new Instance(numberOfVariables_ + 1); //+1 for the objective
			 for(int j=0; j<numberOfVariables_;j++)
				 iExample.setValue((Attribute)solutionParameters.elementAt(j), solution.getDecisionVariables()[j].getValue());
				 
			 //set the objective value
			 iExample.setValue((Attribute)solutionParameters.elementAt(numberOfVariables_), solution.getObjective(i));
			 
			 // add the instance
			 tableOfSolutionsForEachObjective[i].add(iExample);		
		 }
  }
	
	
	@Override
  public void evaluate(SolutionSet solutions) throws JMException, ClassNotFoundException {
		
		for (int i = 0; i < solutions.size(); i++) {
			evaluate(solutions.get(i));
		}
  }

	@Override
  public void evaluate(Solution solution) throws JMException, ClassNotFoundException {
		Object[] output;
		double variance;
		Instance iExample;
		double[] distributionForInstance;
		
		//for every objective...
		for(int i=0; i< numberOfObjectives_; i++)
		{
			iExample = new Instance(numberOfVariables_ + 1); //+1 for the objective
			for(int j=0; j<numberOfVariables_;j++)
				 iExample.setValue((Attribute)solutionParameters.elementAt(j), solution.getDecisionVariables()[j].getValue());
				 
			// Specify that the instance belong to the training set 
			// in order to inherit from the set description
			iExample.setDataset(new Instances("Rel", solutionParameters, capacity));
			
			//approximate
			try {
	      distributionForInstance = tableOfRFClassifiers[i].distributionForInstance(iExample);	      
      } 
			catch (Exception e) { e.printStackTrace(); }
		}
		
	
  }

	private void InitializeTrainingSets() {
		
		solutionParameters = new FastVector(numberOfVariables_ + 1);	//+1 for the objective
    Attribute attribute;
    
   
    //Variables
    for(int i=0; i< numberOfVariables_; i++)
    {
    	attribute = new Attribute(i+ " - input variable");
    	solutionParameters.addElement(attribute);
    }
    
    //objective		   
    attribute = new Attribute("objective");
    solutionParameters.addElement(attribute);
    
    
    // Create an empty training set
    Instances trainingSetOfSolutions = new Instances("Rel", solutionParameters, capacity);         
    // Set class index
    trainingSetOfSolutions.setClassIndex(numberOfVariables_);
	
		for(int i=0;i<numberOfObjectives_; i++)
			tableOfSolutionsForEachObjective[i] = new Instances(trainingSetOfSolutions);

}
	
	private void initializeRandomForests()
	{		
		for(int i=0; i<numberOfObjectives_; i++)
		{
			tableOfRFClassifiers[i] = new RandomForest();
			tableOfRFClassifiers[i].setNumTrees(10);
			tableOfRFClassifiers[i].setSeed((int)Params.getSeed());
		}
	}
} // RF
