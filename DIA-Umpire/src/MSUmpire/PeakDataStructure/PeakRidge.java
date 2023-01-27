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
package MSUmpire.PeakDataStructure;

import java.io.Serializable;

/**
 * Peak ridge class (described in 
 * Tautenhahn, R., Bottcher, C. & Neumann, S. 
 * Highly sensitive feature detection for high resolution LC/MS. 
 * BMC Bioinformatics 9, 504 (2008).
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PeakRidge implements Comparable<PeakRidge>, Serializable {

    public float RT;
    public int lowScale;
    public int ContinuousLevel = 0;
    public float intensity;

    @Override
    public int compareTo(PeakRidge o) {
        return Float.compare(o.RT, RT);
    }
}
