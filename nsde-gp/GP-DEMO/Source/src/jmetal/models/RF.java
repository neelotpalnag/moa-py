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

import java.util.ArrayList;
import java.util.List;

import jmetal.core.*;
import jmetal.init.Params;
import jmetal.util.JMException;

import org.rosuda.JRI.RMainLoopCallbacks;
import org.rosuda.JRI.Rengine;
import org.rosuda.JRI.REXP;

//import hr.irb.fastRandomForest;
/**
 * Class representing model Random Forest
 */
public class RF extends Model {  
  
	int numberOfVariables_;
	int numberOfObjectives_;
	List<Solution> trainingSolutions;
	Rengine rengine;
	int numberOfElementsInATree;
	int noOfTrees;
	
	double[] tableOfSolutionsForEachObjective;
	
	//int capacity = 100;
	
  /** 
   * Constructor.
   * Creates a default instance of the GP model.
   * @throws Exception 
   */
  public RF(int noInput, int noOutput, int nTrees) throws Exception {
  	super(noInput, noOutput);
  	  	
  	numberOfVariables_   = noInput;		
	  numberOfObjectives_  = noOutput;
	  	  
	  trainingSolutions =new ArrayList<Solution>();
	  
	  String[] Rags = {"--vanilla"};
	  rengine = new Rengine(Rags, false, null);		//,,new CallbackListener()
	  rengine.eval("library(randomForest)");
	  rengine.eval(String.valueOf(Params.getSeed()));	  
	  
	  numberOfElementsInATree = 5;
	  
	  noOfTrees = nTrees;
	  
	  
	  //initializeRandomForests();
	  //InitializeTrainingSets();
  } // RF

  @Override
  public void evaluate(SolutionSet solutions) throws JMException, ClassNotFoundException {
		
		for (int i = 0; i < solutions.size(); i++) {
			evaluate(solutions.get(i));
		}
  }

	@Override
  public void evaluate(Solution solution) throws JMException, ClassNotFoundException {
		
		double[] ldObjectiveValuePredictions = new double[numberOfObjectives_];
		double[] ldObjectiveValueStdDev = new double[numberOfObjectives_];
		
		/*
		String rFormulaToInitValues = "c(";
	  //za vse testne podatke zgradimo data frame v R-ju
		for(int i=0;i<numberOfVariables_;i++)
		{
			rFormulaToInitValues = "variablePredict" + i + " <- c(" +solution.getDecisionVariables()[i].getValue() + ")";
			rengine.eval(rFormulaToInitValues);
		}
			
		
		for(int i=0; i< numberOfObjectives_;i++)
		{
			rFormulaToInitValues = "objectivePredict" + i + " <- c(" +solution.getObjective(i) + ")";
			rengine.eval(rFormulaToInitValues);
		}
		
		//create data frames for prediction
		String lsDataFrameVariables = " <- data.frame(";
		//spremenljivke
		for(int i=0; i< numberOfVariables_;i++)		
			lsDataFrameVariables = lsDataFrameVariables + "variablePredict" + i + ", ";		//the name "objective1 was defined earlier		
		  	
	  
		
		//kriteriji
		String lsDataFrame;
		for(int i=0; i< numberOfObjectives_;i++)	
		{
			//for every objective we create seperate data frame (df0, df1,...)
			lsDataFrame = "dfPredict" + i + lsDataFrameVariables + "objectivePredict" + i + ")";
			rengine.eval(lsDataFrame);
		}
		*/
		
		
		
		String lsVariableNames = "";
		String lsDataFrame;
		//spremenljivke
		for(int i=0; i< numberOfVariables_-1;i++)		
			lsVariableNames = lsVariableNames + "variable" + i + "= " + solution.getDecisionVariables()[i].getValue() + ", ";		//the name "objective1 was defined earlier		
		lsVariableNames = lsVariableNames + "variable" + (numberOfVariables_-1) + "= " + solution.getDecisionVariables()[numberOfVariables_-1].getValue();	
		
		for(int i=0; i< numberOfObjectives_;i++)	
		{
			lsDataFrame = "dfm" + i + "<- data.frame(" + lsVariableNames + ", objective" + i + "= " + solution.getObjective(i) + ")";
			rengine.eval(lsDataFrame);
		}
		
		
		
		
		//make predictions
		String lsMakePrediction;
		for(int i=0; i< numberOfObjectives_;i++)
		{
			 lsMakePrediction ="predictions" + i + " <- predict(model" + i + ", newdata=dfm" + i + ", predict.all=TRUE)";	//dfPredict
			 rengine.eval(lsMakePrediction);
		}
		
		//calculate mean 
		String lsgetMean;
		REXP result;
		for(int i=0; i< numberOfObjectives_;i++)
		{
			 lsgetMean ="mean" + i + " <- apply(predictions" + i + "$individual, MARGIN=1, mean)";
			 result = rengine.eval(lsgetMean);
			 ldObjectiveValuePredictions[i] = result.asDoubleArray()[0];   // na niè, ker imamo le eno predikcijo rešitve
			 solution.setObjective(i, ldObjectiveValuePredictions[i]);
		}
		
		//calculate standard deviance 
		String lsgetStdDev;
		for(int i=0; i< numberOfObjectives_;i++)
		{
			lsgetStdDev ="sd" + i + " <- apply(predictions" + i + "$individual, MARGIN=1, sd)";
			result = rengine.eval(lsgetStdDev);
			ldObjectiveValueStdDev[i] = result.asDoubleArray()[0] *2;
			solution.setStandardDeviance(i, ldObjectiveValueStdDev[i]);
		}
		
		
		
	}
	
	
	@Override
  public void update(SolutionSet solutions) {
		
		for(int i=0; i<solutions.size();i++)
			trainingSolutions.add(solutions.get(i));
      
		buildClassifiers();
	
  }

	
	@Override
  public void update(Solution solution) {
		
		trainingSolutions.add(solution);
		
		buildClassifiers();
  }

	@Override
	public void closeRConnection()
	{
		rengine.end();			
	}
	
	private void buildClassifiers() {

		String rFormula;
	  //za vse testne podatke zgradimo data frame v R-ju
		for(int i=0;i<numberOfVariables_;i++)
		{
			rFormula = getColumnOfTrainingVariableValues(i);
			rengine.eval(rFormula);
		}
		for(int i=0; i< numberOfObjectives_;i++)
		{
			rFormula = getColumnOfTrainingObjectivesValues(i);
			rengine.eval(rFormula);
		}
		
		//iz stolpcev spremenljivk in kriterijev zgradimo data frame
		String lsDataFrameVariables = " = data.frame(";
		//spremenljivke
		for(int i=0; i< numberOfVariables_-1;i++)		
			lsDataFrameVariables = lsDataFrameVariables + "variable" + i + ", ";		//the name "objective1 was defined earlier		
		lsDataFrameVariables = lsDataFrameVariables + "variable" + (numberOfVariables_-1) + ", ";	
	  
		
		//kriteriji
		String lsDataFrame;
		for(int i=0; i< numberOfObjectives_;i++)	
		{
			//for every objective we create seperate data frame (df0, df1,...)
			lsDataFrame = "df" + i + lsDataFrameVariables + "objective" + i + ")";
			rengine.eval(lsDataFrame);
		}
			//rFormula = rFormula + "objective" + i + ", ";		//the name "objective1 was defined earlier		
		//rFormula = rFormula + "objective" + (numberOfObjectives_-1) + ")";	
	  
		//create RF models
		String lsCreateModel;
		for(int i=0; i< numberOfObjectives_; i++)
		{
			lsCreateModel = "model" + i + " <- randomForest(objective" + i + " ~.,data=df" + i + ", ntree = "+noOfTrees + ", nodesize=" + numberOfElementsInATree + ")";
			rengine.eval(lsCreateModel);			
		}
		
		
		
  }

	private String getColumnOfTrainingObjectivesValues(int objective) {
		String termInR;
	  
	  termInR ="c(";
	  
	  for(int i=0; i< trainingSolutions.size()-1; i++)
	  	termInR = termInR + trainingSolutions.get(i).getObjective(objective) + ", ";
	  termInR = "objective"+objective+" <-" + termInR + trainingSolutions.get(trainingSolutions.size()-1).getObjective(objective) + ")";
	  
	  return termInR;
  }

	private String getColumnOfTrainingVariableValues(int variable) {
	  String termInR;
	  
	  termInR ="c(";
	  
	  for(int i=0; i< trainingSolutions.size()-1; i++)
	  	termInR = termInR + trainingSolutions.get(i).getDecisionVariables()[variable] + ", ";
	  termInR = "variable"+variable+" <-" + termInR + trainingSolutions.get(trainingSolutions.size()-1).getDecisionVariables()[variable] + ")";
	  
	  return termInR;
  }


	
private static class CallbackListener implements RMainLoopCallbacks {
		
		public void rWriteConsole(Rengine arg0, String arg1, int ii) {
			System.out.print(arg1);
			}
			public void rBusy(Rengine arg0, int arg1) { }
			public String rReadConsole(Rengine arg0, String arg1, int arg2) { return null; }
			public void rShowMessage(Rengine arg0, String arg1) { }
			public String rChooseFile(Rengine arg0, int arg1) { return null; }
			public void rFlushConsole(Rengine arg0) { }
			public void rSaveHistory(Rengine arg0, String arg1) { }
			public void rLoadHistory(Rengine arg0, String arg1) { }
	}
	
	
	/*private void buildClassifiers() throws Exception {

		for(int i=0; i< numberOfObjectives_; i++)
			tableOfRFClassifiersRengine[i].buildClassifier(tableOfSolutionsForEachObjective[i]);		//Classifier cModel = (Classifier)new NaiveBayes();
			  
  }*/
/*
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
	*/
	

	
	/*
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
	      distributionForInstance = tableOfRFClassifiersRengine[i].distributionForInstance(iExample);	      
      } 
			catch (Exception e) { e.printStackTrace(); }
		}
		
	
  }*/
/*
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
			tableOfRFClassifiersRengine[i] = new RandomForest();
			tableOfRFClassifiersRengine[i].setNumTrees(10);
			tableOfRFClassifiersRengine[i].setSeed((int)Params.getSeed());
		}
	}*/
} // RF
