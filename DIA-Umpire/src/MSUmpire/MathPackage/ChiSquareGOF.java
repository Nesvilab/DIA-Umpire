/*
 * This file is part of DIA-Umpire.
 *
 * DIA-Umpire is free software: you can redistribute it and/or modify it under 
 * the terms of the GNU General Public License as published by the Free Software 
 * Foundation, either version 3 of the License, or (at your option) any later 
 * version.
 *
 * DIA-Umpire is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
 * more details.
 *
 * You should have received a copy of the GNU General Public License along 
 * with DIA-Umpire. If not, see <https://www.gnu.org/licenses/>.
 *
 */
package MSUmpire.MathPackage;

import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;
import jsc.distributions.ChiSquared;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class ChiSquareGOF {

    private static ChiSquareGOF models = null;
    public static ChiSquared[] chimodels;
    public static ReadWriteLock lock = new ReentrantReadWriteLock();    
    
    private ChiSquareGOF(int maxpeak) {      
        chimodels = new ChiSquared[maxpeak-1];
        for (int i = 1; i < maxpeak; i++) {
            chimodels[i-1] = new ChiSquared(i);
        }
    }

    public static ChiSquareGOF GetInstance(int maxpeak) {
        if (models == null) {
            lock.writeLock().lock();
            try {
                if (models == null) {
                    models = new ChiSquareGOF(maxpeak);
                }
            } finally {
                lock.writeLock().unlock();
            }
        }
        return models;
    }

    public float GetGoodNessOfFitProb(float[] expected, float[] observed) {
        float gof = 0f;
        int nopeaks = 0;
        for (int i = 0; i < Math.min(observed.length, expected.length); i++) {
            if (observed[i] > 0) {
                float error = expected[i] - observed[i];
                gof += (error * error) / (expected[i] * expected[i]);
                nopeaks++;
            }
        }
        if (Float.isNaN(gof) || nopeaks < 2) {
            return 0f;
        }
        
        if(chimodels[nopeaks-2]==null){
            System.out.println("");
        }
        
        float prob = 1 - (float) chimodels[nopeaks - 2].cdf(gof);
        return prob;
    }
}