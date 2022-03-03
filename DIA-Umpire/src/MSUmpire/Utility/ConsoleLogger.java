/* 
 * Author: Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 *             Nesvizhskii Lab, Department of Computational Medicine and Bioinformatics, 
 *             University of Michigan, Ann Arbor
 *
 * Copyright 2014 University of Michigan, Ann Arbor, MI
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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
