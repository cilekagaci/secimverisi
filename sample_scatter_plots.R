library('tidyverse')


data = read_csv('genel25_26_referandum2017_27meclis_cb2018.csv') 

# 2015 Kasim - 2018 Haziran MHP oy oranlari degisimi
data2018 %>%
  ggplot(aes(x = MHP.Kasim / total.valid.vote.count.Kasim, 
             y = MHP.2018meclis / total.valid.vote.count.2018meclis)) +
  geom_point()


# 2015 Kasim - 2018 Haziran gecersiz oy oranlari degisimi
data2018 %>%
  ggplot(aes(x = invalid.vote.count.Kasim / voter.count.Kasim, 
             y = invalid.vote.count.2018meclis / voter.count.2018meclis)) +
  geom_point()


# 2018 Haziran Cumhur Ittifaki - Erdogan oylari karsilastirmasi
data2018 %>%
  ggplot(aes(x = (AKP.2018meclis + MHP.2018meclis) / total.valid.vote.count.2018meclis, 
             y = ERDOGAN.2018cb / total.valid.vote.count.2018cb)) +
  geom_point()


# Erdogan'in Cumhur ittifakinin en az %5 ustunde oy aldigi ilcelerde fark ve Bagimsiz+HUDAPAR+SP oy oranlari
data2018 %>%
  mutate(fark = ERDOGAN.2018cb / total.valid.vote.count.2018cb -
           (AKP.2018meclis + MHP.2018meclis) / total.valid.vote.count.2018meclis) %>%
  filter(fark >= 0.03) %>%
  ggplot(aes(x = (SP.2018meclis + HUDAPAR.2018meclis +  BGMSZ.2018meclis) / total.valid.vote.count.2018meclis, 
             y = fark,
             size = voter.count.2018meclis)) +
  geom_point() +
  geom_abline(linetype = 'dashed')
