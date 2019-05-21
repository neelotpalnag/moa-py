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
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import com.mathworks.toolbox.javabuilder.MWCellArray;
import com.mathworks.toolbox.javabuilder.MWException;
import com.mathworks.toolbox.javabuilder.MWNumericArray;

import sun.misc.BASE64Encoder;

import jmetal.core.*;
import jmetal.util.JMException;
import gpml.GPML;

/**
 * Class representing model GP (Gaussian Process)
 */
public class GP extends Model {  
  
	private GPML gp_;
	private boolean fixed_;
	private boolean reset_;
	private int hypLength_;
	private int activeSetLimit_;
	private MWCellArray covfunc_;
	private MWNumericArray[] hyp_;
	private List<double[]> inputs_;
	private List<double[]> outputs_;
	private List<String> data_;

  /** 
   * Constructor.
   * Creates a default instance of the GP model.
   * @throws ClassNotFoundException
   * @throws MWException 
   */
  public GP(int noInput, int noOutput) throws Exception {
  	this(noInput, noOutput, -1, false, "covSEard", null);
  }
  
  /** 
   * Constructor.
   * Creates a default instance of the GP model.
   * @throws Exception 
   */
  public GP(int noInput, int noOutput, int activeSetLimit, boolean reset, String covFuncM, double[][] hyp) throws Exception {
  	super(noInput, noOutput);
  	MWCellArray covfunc;
  	
  	reset_ = reset;
  	modelName_ = "GP";
  	// covariance function
  	covfunc = new MWCellArray(2, 1);
  	covfunc.set(new int[]{1,1}, covFuncM);
		covfunc.set(new int[]{2,1}, "covNoise");
		covfunc_ = new MWCellArray(2, 1);
		covfunc_.set(new int[]{1,1}, "covSum");
		covfunc_.set(new int[]{2,1}, covfunc);
  	// GP models
		activeSetLimit_ = activeSetLimit;
  	gp_ = new GPML();
  	gp_.init(); // random seed set
  	hyp_ = new MWNumericArray[noOutput];
  	// get number of hyperparameter according to covariance function and number of inputs
  	hypLength_ = ((MWNumericArray)gp_.hypLen(1, covfunc_, noInput)[0]).getInt();
  	// set hyperparameter values
  	if (hyp == null) {
  		hyp = new double[hypLength_][];
	  	for (int i = 0; i < hypLength_; i++)
			  hyp[i] = new double[]{0};
			for (int i = 0; i < noOutput; i++)
	  		hyp_[i] = new MWNumericArray(hyp); // default hyperparameter values
			fixed_ = false;
  	}
  	else {
  		if (hyp.length != hypLength_ || hyp[0].length != noOutput)
  			throw new Exception("Dimenzija podanih hiperparametrov ni OK");
  		double[][] hypt;
			for (int i = 0; i < noOutput; i++) {
				hypt = new double[hypLength_][];
				for (int j = 0; j < hypLength_; j++)
					hypt[j] = new double[]{hyp[j][i]};
	  		hyp_[i] = new MWNumericArray(hypt); // default hyperparameter values
			}
  		fixed_ = true;
  	}
  	// set data
  	inputs_ = new ArrayList<double[]>();
  	outputs_ = new ArrayList<double[]>();
  	data_ = new ArrayList<String>();
  } // GP

	@Override
  public void update(SolutionSet solutions) throws JMException, ClassNotFoundException {
		Solution solution;
		
		for (int i = 0; i < solutions.size(); i++) {
			solution = solutions.get(i);
			append(solution);
		}
		try {
			if (!fixed_)
				train();
				//trainParallel();
    }
		catch (Exception e) {
	    e.printStackTrace();
    }
  }

	@Override
  public void update(Solution solution) throws JMException, ClassNotFoundException {
		append(solution);
		try {
			if (!fixed_)
				train();
				//trainParallel();
    }
		catch (Exception e) {
	    e.printStackTrace();
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
		
		for (int i = 0; i < numberOfOutputs_; i++) {
			try {
	      output = gp_.gpr(2, hyp_[i], covfunc_, toArray(inputs_), toArray(outputs_, i), toArray(solution.getDecisionVariables()));
	      if (solution.numberOfFeatures() > 0)
	      	solution.setFeature(i, ((MWNumericArray)output[0]).getDouble());
	      else
	      	solution.setObjective(i, ((MWNumericArray)output[0]).getDouble());
	      variance = ((MWNumericArray)output[1]).getDouble();
	      solution.setStandardDeviance(i, Math.sqrt(variance) * 2);
	      solution.setExactllyEvaluated(false);
      }
			catch (MWException e) {
	      e.printStackTrace();
      }
		}
  }

	private void append(Solution solution) throws JMException {
		double[] objectives;
		double[] variables;
		
		// variables
		variables = new double[solution.numberOfVariables()];
		for (int i = 0; i < solution.numberOfVariables(); i++)
			variables[i] = solution.getDecisionVariables()[i].getValue();
		// objectives
		if (solution.numberOfFeatures() > 0) {
			objectives = new double[solution.numberOfFeatures()];
			for (int i = 0; i < objectives.length; i++)
				objectives[i] = solution.getFeature(i);
		}
		else {
			objectives = new double[solution.numberOfObjectives()];
			for (int i = 0; i < objectives.length; i++)
				objectives[i] = solution.getObjective(i);
		}
		// check if already exists
		if (data_.contains(toHash(variables, objectives)))
			return;
		// append
		inputs_.add(variables);
		outputs_.add(objectives);
		data_.add(toHash(variables, objectives));
		if (activeSetLimit_ > 0 && inputs_.size() > activeSetLimit_) {
			inputs_.remove(0);
			outputs_.remove(0);
			data_.remove(0);
		}
	}

	private String toHash(double[] variables, double[] objectives) {
		BASE64Encoder base64;
		StringBuilder concat;
		MessageDigest md;
		
		concat = new StringBuilder();
		base64 = new BASE64Encoder();
		try {
			md = MessageDigest.getInstance("SHA-1");
		}
		catch (NoSuchAlgorithmException e) {
			// tale more obstajat
			md = null;
		}
		for (int i = 0; i < variables.length; i++)
			concat.append(Long.toHexString(Double.doubleToLongBits(variables[i])) + "|");
		concat.append("-|");
		for (int i = 0; i < objectives.length; i++)
			concat.append(Long.toHexString(Double.doubleToLongBits(objectives[i])) + "|");
		return base64.encode(md.digest(concat.toString().getBytes()));
	}
	
	private void train() throws MWException {
		Object[] output;
		double[][] hyp;
		
		for (int i = 0; i < numberOfOutputs_; i++) {
			if (reset_) {
				//output = gp_.train(2, covfunc_, toArray(inputs_), toArray(outputs_,i));
	  		hyp = new double[hypLength_][];
		  	for (int j = 0; j < hypLength_; j++)
				  hyp[j] = new double[]{0};
				output = gp_.train(2, covfunc_, toArray(inputs_), toArray(outputs_,i), new MWNumericArray(hyp));
			}
			else
				output = gp_.train(2, covfunc_, toArray(inputs_), toArray(outputs_,i), hyp_[i]);
			hyp_[i] = (MWNumericArray)output[0];
		}
	}	

	private void trainParallel() throws MWException, ExecutionException, InterruptedException {
		List<Callable<MWNumericArray>> tasks;
		List<Future<MWNumericArray>> futures;
		ExecutorService exec;
		
		tasks = new ArrayList<Callable<MWNumericArray>>();
		exec = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
		for (int i = 0; i < numberOfOutputs_; i++) {
			tasks.add(new GPTrainCallable(gp_, covfunc_, toArray(inputs_), toArray(outputs_, i)));
		}
		
		futures = exec.invokeAll(tasks);
		
		for (int i = 0; i < futures.size(); i++)
      hyp_[i] = futures.get(i).get();
	}
	
	private double[][] toArray(List<double[]> list) {
		double[][] array;
		
		array = new double[list.size()][];
		for (int i = 0; i < array.length; i++)
			array[i] = list.get(i);
		return array;
	}

	private double[][] toArray(List<double[]> list, int index) {
		double[][] array;
		
		array = new double[list.size()][];
		for (int i = 0; i < array.length; i++)
			array[i] = new double[]{list.get(i)[index]};
		return array;
	}

	private double[] toArray(Variable[] variables) throws JMException {
		double[] array;
		
		array = new double[variables.length];
		for (int i = 0; i < array.length; i++)
			array[i] = variables[i].getValue();
		return array;
	}
	
	private String toString(List<double[]> list) {
		return toString(list, 0);
	}
	
	private String toString(List<double[]> list, int start) {
		StringBuilder sb;
		double[] values;
		
		sb = new StringBuilder();
		for (int i = start; i < list.size(); i++) {
			values = list.get(i);
			for (int j = 0; j < values.length; j++) {
				sb.append(values[j]);
				sb.append("\t");
			}
			sb.append("\n");
		}
		return sb.toString();
	}
	
} // GP

class GPTrainCallable implements Callable<MWNumericArray> {
	
	private GPML gp_;
	private MWCellArray covfunc_;
	private double[][] inputs_;
	private double[][] outputs_;
	
	public GPTrainCallable(GPML gp, MWCellArray covfunc, double[][] inputs, double[][] outputs) {
		gp_ = gp;
		covfunc_ = covfunc;
		inputs_ = inputs;
		outputs_ = outputs;
	}
	
	@Override
	public MWNumericArray call() throws MWException {
		Object[] output;
		
		System.out.println("Started");
		output = gp_.train(2, covfunc_, inputs_, outputs_);
		System.out.println("Finished");
		return (MWNumericArray)output[0];
	}
	
}