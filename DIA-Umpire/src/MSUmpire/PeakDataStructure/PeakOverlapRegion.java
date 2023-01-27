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
 * Coeluting peak pair
 * @author Chih-Chiang Tsou
 */
public class PeakOverlapRegion implements Serializable{
    private static final long serialVersionUID = 2654984654657411L;

    public int PeakCurveIndexA;
    public int PeakCurveIndexB;
    public float Correlation;
}
