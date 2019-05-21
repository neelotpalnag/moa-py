//  NN.java
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

import weka.classifiers.functions.MultilayerPerceptron;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import jmetal.core.*;
import jmetal.init.Params;
import jmetal.util.JMException;
import jmetal.util.Configuration.*;

/**
 * Class representing model NN (Neural Network)
 */
public class NN extends Model {  
    
	private int noOfInput;
	private int noOfOutput;
	MultilayerPerceptron [] tableOfNeuralNetworks;
	Instances [] iExamples;
  /** 
   * Constructor.
   * Creates a default instance of the NN model.
   */
  public NN(int noInput, int noOutput) throws ClassNotFoundException {
  	super(noInput, noOutput);
  	noOfInput = noInput;
  	noOfOutput = noOutput;
  	
  	modelName_ = "NN";
  	
  	tableOfNeuralNetworks = new MultilayerPerceptron[noOutput];
  	iExamples = new Instances[noOutput];
  	for(int i=0; i< noOutput; i++)
  	{
  		tableOfNeuralNetworks[i] = new MultilayerPerceptron();
  		iExamples[i] = getInitializedInstances(noInput);
  	}
  	
		double evaluatedValue;
		double ldDistanceToFisibility;
  	
  	
  	
  	
  } // NN
      
  /** 
  * Updates a model
  * @param input The input...
  * @param target The target... 
  * @throws JMException 
  */
  public void update(double[] input, double[] target) throws JMException {
  } // update

  /** 
  * Make a prediction
  * @param input The input...
  * @throws JMException 
  */
  public double[] predict(double[] input) throws JMException {
	  return null;
  } // predict

	@Override
  public void update(SolutionSet solutions) throws JMException, ClassNotFoundException {
  }

	@Override
  public void update(Solution solution) throws JMException, ClassNotFoundException {
  }

	@Override
  public void evaluate(SolutionSet solutions) throws JMException, ClassNotFoundException {
  }

	@Override
  public void evaluate(Solution solution) throws JMException, ClassNotFoundException {
		
		Instance iExample = new Instance(Params.getNumberOfVariables() + 1);
		
		for(int k=0; k<noOfOutput;k++)      				
			iExample.setValue(k, solution.getDecisionVariables()[k].getValue());
		
  }

	private Instances getInitializedInstances(int noInput) {

    FastVector attr = new FastVector(noInput + 1);	//+1 for the criterion
    Attribute Attribute1;
    for(int i=0; i< noInput; i++)
    {
    	Attribute1 = new Attribute(i+ "- input variable");
    	attr.addElement(Attribute1);
    }
    
    Attribute1 = new Attribute("- class");
	attr.addElement(Attribute1);
	
	

Instances iExamples = new Instances("Examples", attr, 0); 
iExamples.setClassIndex(Params.getNumberOfVariables());
return iExamples;
}
	
} // NN