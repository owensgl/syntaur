
data <- read_anchorwave("Hneg2702r1.0_H1_newname_vs_Ha412_anchorwave.tab_tab_fp.txt")

data <- read_anchorwave("data/sunflower.artichoke.anchorwave.blocks.txt")
matches <- find_matches(data, 10000000)

test <- data[which(data$qid == "Chr16" & data$rid == "Chr17"),c(1:4)]
test <- data[which(data$qid == "CM060405.1" & data$rid == "Ha412HOChr10"),c(1:4)]

test_all_flips(test)
basic_alignment_plot(test)

result <- multi_iterate(test,iterations=100,verbose=T)
list_flipped_windows(test, result)

plotting_result <- plot_flip_schedule(test, result)

plot_final_series(plotting_result, result)

find_all_inversions("Hneg2702r1.0_H1_newname_vs_Ha412_anchorwave.tab_tab_fp.txt", "Hneg2702r1.0_H1",iterations=100)
find_all_inversions("Hneg2702r1.0_H2_newname_vs_Ha412_anchorwave.tab_tab_fp.txt", "Hneg2702r1.0_H2",iterations=100)

result <- test_all_flips(test)

test %>%
  ggplot(.,aes(x=rs,y=qs)) + geom_point()


test_sim <- simulate_inversions(inversions=2)
multi_iterate(test_sim[[2]],iterations=100)

simulate_inversions(inversions=5) %>%
  ggplot(.) +
  geom_segment(ggplot2::aes(x = rs, xend = re, y = qs, yend = qe, color = rs),size=2) +
  xlab("Reference") + ylab("Query") +
  scale_color_viridis_c(name="Reference position") +
  ggtitle("Simulated inversions: 5")

