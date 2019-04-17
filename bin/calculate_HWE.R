file <- commandArgs()[4:4]
table <- read.table(file, row.names=1, header=TRUE)
row.fisher <- function(row, alt='two.sided', cnf=0.95, hybrid=TRUE) {
f <- fisher.test(matrix(row, nrow=3), alternative=alt, conf.level=cnf, hybrid=TRUE)
return(c(row, p_val=f$p.value, or=f$estimate[[1]], or_ll=f$conf.int[1], or_ul=f$conf.int[2]))
}
p <- t(apply(table, 1, row.fisher))
write.table(p, "IMSindel_deletions_HWE_results.txt", sep="\t")
