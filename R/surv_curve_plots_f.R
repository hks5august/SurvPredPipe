#' This function generate individual survival curve for patients. Moreover, user can also highlight one specific sample by providing sample IDs
#' @param Surv_curve_data :args1 - Survival probabilty data for survival curve of patients 
#' @param selected_sample:args2 - ID of a specific sample user wants to highlights on the plot against the background of curves of other samples 
#' @examples
#' SurvPredPipe::surv_curve_plots_f(Surv_curve_data="survCurves_data.txt", selected_sample="TCGA-DH-5140-01")
#' usgae: surv_curve_plots_f(Surv_curve_data, selected_sample)
#' @export

surv_curve_plots_f <- function(Surv_curve_data, selected_sample )   {


#survCurves_data <- read.table("survCurves2_data.txt", header = TRUE, sep = "\t",  check.names = FALSE)
survCurves_data <- read.table(Surv_curve_data, header = TRUE, sep = "\t",  check.names = FALSE)


#reshape the data in the form of long matrix
survCurves_m <- melt(survCurves_data , id='time_point')

#head(survCurves_m)

colnames(survCurves_m) <- c("Time", "Patient", "Value")
head(survCurves_m)

Surv_curv_plot_all_pat <- ggplot(survCurves_m, aes(x = Time, y=Value, group = Patient, color=Patient)) +
  geom_line() + labs(x="Time in Months",y="Survival Probability")+ 
  ggtitle("Survival Curves for Patients")  +
  geom_hline(yintercept = 0.5, color="black", linetype="dashed") + #add dashed line corresonds to 0.5 probability
  theme( legend.position="bottom", legend.box = "vertical", legend.title = element_text(), legend.key = element_blank(),
         legend.key.size = unit(3, 'mm'),legend.text = element_text(size=4))  + guides(color = guide_legend(title = "Patients")) 

#Surv_curv_plot_all_pat

jpeg(file="survival_curves_all_patients_with_PI_clin.jpeg", units="in", width=10, height=10, res=300)
print(Surv_curv_plot_all_pat)
dev.off()
 


head(survCurves_m)

#Selected_patient<- as.character("TCGA-DU-A6S3-01")
Selected_patient<- selected_sample

Surv_curv_plot_all_pats_with_highlighting_one_pat <- ggplot(survCurves_m, aes(x = Time, y = Value, linetype = Patient, color = Selected_patient)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  ggtitle("Survival Curves for Patients with Highlightited Patient") +
  geom_hline(yintercept = 0.5, color = "red", linetype = "dashed") +
  
  # Set linetype for all lines based on Patient
  geom_line(size = 0.5) +
  
  # Highlight one patient with a specific color
  scale_linetype_manual(values = rep("solid", length(unique(survCurves_m$Patient)))) +
  scale_color_manual(values = c("black", "yellow")) +
  geom_line(aes(colour = "yellow"), size=0.5,data = ~ subset(survCurves_m, Patient == Selected_patient )) +
  
  labs(x = "Time in Months", y = "Survival Probability") +
  
  theme(
    legend.position = "bottom",
    legend.title = element_text(),
    legend.key = element_blank(),
    legend.key.size = unit(3, 'mm'),
    legend.text = element_text(size = 4)
  ) +
  
  guides(linetype = guide_legend(title = "Patients"), color = guide_legend(title = "Selected Patient"))

Surv_curv_plot_all_pats_with_highlighting_one_pat 

jpeg(file="survival_curves_for_all_patients_with_highlighting_one_pat_PI_clin.jpeg", units="in", width=10, height=10, res=300)
print(Surv_curv_plot_all_pats_with_highlighting_one_pat )
dev.off()



}



