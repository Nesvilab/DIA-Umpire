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
package MSUmpire.Utility;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.core.LoggerContext;
import org.apache.logging.log4j.core.appender.ConsoleAppender;
import org.apache.logging.log4j.core.appender.FileAppender;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.filter.ThresholdFilter;
import org.apache.logging.log4j.core.layout.PatternLayout;
/**
 *
 * @author Chih-Chiang Tsou
 */
public class ConsoleLogger {
    public static void SetConsoleLogger(Level level) {
        final LoggerContext ctx = (LoggerContext) LogManager.getContext(false);
        ctx.getConfiguration().getLoggerConfig(LogManager.ROOT_LOGGER_NAME).removeAppender("Console");
        ctx.getConfiguration().getLoggerConfig(LogManager.ROOT_LOGGER_NAME).setLevel(Level.ALL);
        ConsoleAppender ca = ConsoleAppender.newBuilder()
                .setFilter(ThresholdFilter.createFilter(level, null, null))
                .setName("ConsoleLogger_Info")
                .setLayout(PatternLayout.newBuilder().withPattern("%d %-5p [%t{1}] %m%n").build())
                .build();
        ca.start();
        ctx.getConfiguration().addAppender(ca);
        ctx.getRootLogger().addAppender(ctx.getConfiguration().getAppender(ca.getName()));
        ctx.updateLoggers();
    }

    public static void SetFileLogger(Level level, String filename) {
        final LoggerContext ctx = (LoggerContext)LogManager.getContext(false);
        FileAppender fa = FileAppender.newBuilder()
                .setName("FileLogger_Debug")
                .withFileName(filename)
                .setLayout(PatternLayout.newBuilder().withPattern("%d %-5p [%t{1}] %m%n").build())
                .setFilter(ThresholdFilter.createFilter(level, null, null))
                .withAppend(false)
                .build();
        fa.start();
        ctx.getConfiguration().addAppender(fa);
        ctx.getRootLogger().addAppender(ctx.getConfiguration().getAppender(fa.getName()));
        ctx.updateLoggers();
    }

}
