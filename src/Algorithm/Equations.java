/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Algorithm;

import Jama.Matrix;

/**
 *
 * @author Wai Pai Lee
 */
public interface Equations {
    public double[] f(double[] guess);
    public double valF(double[] guess, int eqnNum);
    public double derivative(double [] guess, int row, int col);
    public Matrix jacobian(double [] guess);
    public Matrix jacobianT(double [] guess);
    public int getEqnNum();
    public int getVariateNum();
}
