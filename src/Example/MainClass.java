/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Example;

import Algorithm.RootSolver;

/**
 *
 * @author Wai Pai Lee
 */
public class MainClass {
    public static void main(String[] args) {
        Example exp = new Example();
        
        RootSolver m = new RootSolver();
        
        double [] guess = {1,0.1,0.1};
        double [] ans;
        
        ans = m.solve(exp, guess, 1);
        
        System.out.println(ans[0]);
        System.out.println(ans[1]);
        System.out.println(ans[2]);
    }
}
