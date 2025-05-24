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
package MSUmpire.BaseDataStructure;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class UmpireInfo {

    private static UmpireInfo umpireInfo = null;
    public String Version = "v2.3.3";

    private UmpireInfo() {

    }

    public static UmpireInfo GetInstance() {
        if (umpireInfo == null) {
            umpireInfo = new UmpireInfo();
        }
        return umpireInfo;
    }

}
