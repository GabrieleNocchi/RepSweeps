library(dplyr)

genes <- read.table("tmp_file.txt", h = TRUE)

ortho_minimum <- genes %>% group_by(orthogroup) %>% slice(which.min(mean_emp_p))

ortho_minimum <- data.frame(ortho_minimum)

gabriele_dunnsidak <- function(x, y){
    1 - ((1-x)^y)
}


dunnsidak_orthopvalues <- gabriele_dunnsidak(ortho_minimum$mean_emp_p,ortho_minimum$ortho_size)


ortho_adjusted <- cbind(ortho_minimum, dunnsidak_orthopvalues)


write.table(ortho_adjusted, file = "final_orthogroup.txt", quote = FALSE, row.names = FALSE,col.names=TRUE, sep = "\t")
