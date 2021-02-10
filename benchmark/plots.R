library(dplyr)
library(ggplot2)

package_colors <- c(
  ape = "#ef798a", Biostrings = "#faa613", tidysq = "#688e26",
  seqinr = "#732c2c", bioseq = "#39a27a", `plik`= "#050413")

results_fs <- results %>%
  group_by(alph_size, sq_len, num_sq) %>%
  slice(1) %>%
  ungroup() %>%
  select(alph_size, sq_len, num_sq, obj_size = file_size) %>%
  mutate(package = "plik") %>%
  mutate(num_sq = factor(num_sq, levels = sort(unique(num_sq)),
                         labels = paste0("Liczba sekwencji:\n", sort(unique(num_sq)))))

melt_res <- results %>% 
  group_by(package, sq_len, num_sq, type, file_size) %>% 
  summarise(obj_size = median(obj_size),
            time_value = median(time_value)) %>% 
  ungroup() %>% 
  mutate(num_sq = factor(num_sq, levels = sort(unique(num_sq)),
                         labels = paste0("Liczba sekwencji:\n", sort(unique(num_sq)))))


filter(melt_res, type == "read") %>% 
  select(sq_len, obj_size, package, num_sq) %>%
  bind_rows(results_fs) %>%
  ggplot(aes(x = sq_len, y = obj_size, color = package, linetype = package == "plik")) +
  #geom_point(aes(y = file_size, shape = "File size"), color = "black", size = 4, shape = 15) +
  geom_point(shape = 0, size = 3) +
  geom_line(size = 1) +
  #geom_line(aes(y = file_size), color = "black", linetype = "dashed") +
  facet_wrap(~num_sq, scales = "free_y") +
  scale_color_manual(values = package_colors, guide = guide_legend(override.aes = list(color = package_colors), )) +
  scale_y_continuous("Rozmiar obiektu w bajtach") +
  scale_x_continuous("Średnia długość sekwencji") +
  ggtitle("Porównanie rozmiarów obiektów i rozmiarów pliku") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom", 
        legend.key.width = unit(1, "cm"),
        plot.title = element_text(hjust = 0.5)) +
  guides(linetype = "none",
         color = guide_legend("Format danych", override.aes = list(linetype = c(1,1,2,1,1,1), shape = 0)))


filter(melt_res, type == "read", package != "seqinr") %>% 
  select(sq_len, obj_size, package, num_sq) %>%
  bind_rows(results_fs) %>%
  ggplot(aes(x = sq_len, y = obj_size, color = package, linetype = (package == "plik"))) +
  geom_point(shape = 0, size = 3) +
  geom_line(size = 0.7) +
  #geom_line(aes(y = file_size), color = "black", linetype = "dashed") +
  facet_wrap(~num_sq, scales = "free_y") +
  scale_color_manual(values = package_colors[-4]) +
  scale_y_continuous("Rozmiar obiektu w bajtach") +
  scale_x_continuous("Średnia długość sekwencji") +
  ggtitle("Porównanie rozmiarów obiektów i rozmiarów pliku (bez seqinr)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom", 
        legend.key.width = unit(1, "cm"),
        plot.title = element_text(hjust = 0.5)) +
  guides(linetype = "none",
         color = guide_legend("Format danych", override.aes = list(linetype = c(1,1,2,1,1), shape = 0)))


filter(melt_res, type == "read") %>% 
  ggplot(aes(x = sq_len, y = time_value, color = package)) +
  geom_point() +
  geom_line() +
  facet_wrap( ~ num_sq, scales = "free_y") +
  scale_color_manual(values = package_colors) +
  scale_x_continuous("Średnia długość sekwencji") +
  scale_y_continuous("Czas odczytania") +
  ggtitle("Czas odczytywania z pliku") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom", 
        legend.key.width = unit(1, "cm"),
        plot.title = element_text(hjust = 0.5))


filter(melt_res, type == "char") %>% 
  ggplot(aes(x = sq_len, y = time_value, color = package)) +
  geom_point() +
  geom_line() +
  facet_wrap( ~ num_sq, labeller = label_both, scales = "free_y") +
  scale_color_manual(values = package_colors) +
  scale_x_continuous("Średnia długość sekwencji") +
  scale_y_continuous("Czas konwersji") +
  ggtitle("Czas konwertowania na ciąg znaków") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom", 
        legend.key.width = unit(1, "cm"),
        plot.title = element_text(hjust = 0.5))


filter(melt_res, type == "cons") %>% 
  ggplot(aes(x = sq_len, y = time_value, color = package)) +
  geom_point() +
  geom_line() +
  facet_wrap( ~ num_sq, labeller = label_both, scales = "free_y") +
  scale_color_manual(values = package_colors) +
  scale_x_continuous("Średnia długość sekwencji") +
  scale_y_continuous("Czas konstrukcji") +
  ggtitle("Czas konstrukcji z ciągu znaków") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom", 
        legend.key.width = unit(1, "cm"),
        plot.title = element_text(hjust = 0.5))


filter(melt_res, type == "tran") %>% 
  ggplot(aes(x = sq_len, y = time_value, color = package)) +
  geom_point() +
  geom_line() +
  facet_wrap( ~ num_sq, labeller = label_both, scales = "free_y") +
  scale_color_manual(values = package_colors) +
  scale_x_continuous("Średnia długość sekwencji") +
  scale_y_continuous("Czas translacji") +
  ggtitle("Czas translacji kodonów") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom", 
        legend.key.width = unit(1, "cm"),
        plot.title = element_text(hjust = 0.5))

filter(melt_res, type == "tran", !package %in% c("seqinr", "bioseq")) %>% 
  ggplot(aes(x = sq_len, y = time_value, color = package)) +
  geom_point() +
  geom_line() +
  facet_wrap( ~ num_sq, labeller = label_both, scales = "free_y") +
  scale_color_manual(values = package_colors) +
  scale_x_continuous("Średnia długość sekwencji") +
  scale_y_continuous("Czas translacji") +
  ggtitle("Czas translacji kodonów (bez bioseq i seqinr)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom", 
        legend.key.width = unit(1, "cm"),
        plot.title = element_text(hjust = 0.5))
