context("operation on bytes, packing, unpacking and recoding")

test_that("integers are packed correctly with alph_size = 2, in_len <= 4 ", {
  expect_equal(pack(c(1), 2), as.raw(64))
  expect_equal(pack(c(1, 1), 2), as.raw(80))
  expect_equal(pack(c(1, 1, 1), 2), as.raw(84))
  expect_equal(pack(c(1, 1, 1, 1), 2), as.raw(85))
  
  expect_equal(pack(c(2), 2), as.raw(128))
  expect_equal(pack(c(2, 1), 2), as.raw(144))
  expect_equal(pack(c(2, 2), 2), as.raw(160))
  expect_equal(pack(c(2, 2, 1), 2), as.raw(164))
  expect_equal(pack(c(1, 1, 2), 2), as.raw(88))
  expect_equal(pack(c(2, 1, 2, 1), 2), as.raw(153))
  
  expect_equal(pack(c(3), 2), as.raw(192))
  expect_equal(pack(c(3, 3), 2), as.raw(240))
  expect_equal(pack(c(3, 3, 3, 3), 2), as.raw(255))
  expect_equal(pack(c(1, 2, 2, 3), 2), as.raw(107))
  expect_equal(pack(c(2, 3, 3, 1), 2), as.raw(189))
})

test_that("integers are packed correctly with alph_size = 2, in_len > 4 ", {
  expect_equal(pack(c(1, 1, 1, 1, 1, 1), 2), as.raw(c(85, 80)))
  expect_equal(pack(c(1, 1, 1, 1, 1, 1, 1, 1), 2), as.raw(c(85, 85)))
  expect_equal(pack(c(2, 2, 2, 2, 2, 2, 2), 2), as.raw(c(170, 168)))
  expect_equal(pack(rep(3, 16), 2), as.raw(rep(255, 4)))
  expect_equal(pack(c(1, 1, 3, 2, 3, 1, 3, 2, 3, 1, 3, 2, 3, 1), 2), 
               as.raw(c(94, 222, 222, 208)))
  expect_equal(pack(c(3, 2, 1, 2, 3, 3, 2, 3, 3, 3, 2, 1, 2, 3, 2, 3, 1, 1, 3, 1, 1), 2),
               as.raw(c(230, 251, 249, 187, 93, 64)))
  expect_equal(pack(c(1, 1, 2, 3, 3, 1, 3, 1, 3, 3, 2, 3, 1, 2, 3, 3, 3, 3, 2, 1, 2, 3, 2, 2), 2),
               as.raw(c(91, 221, 251, 111, 249, 186)))
})

test_that("integers are packed correctly with alph_size = 3, in_len <= 2 ", {
  expect_equal(pack(c(1), 3), as.raw(32))
  expect_equal(pack(c(3), 3), as.raw(96))
  expect_equal(pack(c(6), 3), as.raw(192))
  expect_equal(pack(c(2), 3), as.raw(64))
  expect_equal(pack(c(7), 3), as.raw(224))
  expect_equal(pack(c(6, 3), 3), as.raw(204))
  expect_equal(pack(c(7, 2), 3), as.raw(232))
  expect_equal(pack(c(1, 4), 3), as.raw(48))
  expect_equal(pack(c(3, 4), 3), as.raw(112))
  expect_equal(pack(c(7, 7), 3), as.raw(252))
})

test_that("integers are packed correctly with alph_size = 3, in_len > 2 ", {
  expect_equal(pack(c(1, 1, 1, 1, 1, 1, 1), 3), as.raw(c(36, 146, 72)))
  expect_equal(pack(rep(7, 8), 3), as.raw(c(255, 255, 255)))

  expect_equal(pack(c(5, 6, 1), 3), as.raw(c(184, 128)))
  expect_equal(pack(c(5, 3, 6), 3), as.raw(c(175, 0)))
  expect_equal(pack(c(3, 7, 7), 3), as.raw(c(127, 128)))
  expect_equal(pack(c(6, 7, 3), 3), as.raw(c(221, 128)))
  expect_equal(pack(c(3, 1, 1), 3), as.raw(c(100, 128)))
  expect_equal(pack(c(2, 1, 2, 4, 6, 1, 4, 3), 3), as.raw(c(69, 76, 99)))
  expect_equal(pack(c(1, 2, 1, 6, 6, 3, 3, 1), 3), as.raw(c(40, 236, 217)))
  expect_equal(pack(c(3, 6, 3, 3, 3, 4, 5, 2), 3), as.raw(c(121, 183, 42)))
  expect_equal(pack(c(5, 4, 3, 4, 1, 3, 2, 1, 6, 5, 3, 7, 1, 4), 3), as.raw(c(177, 194, 209, 213, 243, 0)))
  expect_equal(pack(c(2, 2, 1, 3, 5, 6, 6, 7, 6, 4, 2, 6, 5, 4, 7, 4, 1, 7, 2, 1), 3), as.raw(c(72, 187, 183, 209, 107, 60, 61, 16)))
  expect_equal(pack(c(4, 4, 1, 5, 6, 4, 6, 3, 1, 2, 6, 2, 4), 3), as.raw(c(144, 221, 51, 43, 40)))
})

test_that("integers are packed correctly with alph_size = 4, in_len <= 2", {
  expect_equal(pack(c(1), 4), as.raw(16))
  expect_equal(pack(c(13), 4), as.raw(208))
  expect_equal(pack(c(3), 4), as.raw(48))
  expect_equal(pack(c(15), 4), as.raw(240))
  expect_equal(pack(c(12), 4), as.raw(192))
  expect_equal(pack(c(10, 1), 4), as.raw(161))
  expect_equal(pack(c(4, 1), 4), as.raw(65))
  expect_equal(pack(c(5, 5), 4), as.raw(85))
  expect_equal(pack(c(10, 12), 4), as.raw(172))
  expect_equal(pack(c(6, 10), 4), as.raw(106))
  expect_equal(pack(c(13, 10), 4), as.raw(218))
  expect_equal(pack(c(4, 6), 4), as.raw(70))
  expect_equal(pack(c(1, 12), 4), as.raw(28))
  expect_equal(pack(c(5, 14), 4), as.raw(94))
  expect_equal(pack(c(4, 14), 4), as.raw(78))
  expect_equal(pack(c(3, 10), 4), as.raw(58))
  expect_equal(pack(c(12, 6), 4), as.raw(198))
  expect_equal(pack(c(15, 15), 4), as.raw(255))
}) 

test_that("integers are packed correctly with alph_size = 5, in_len > 2", {
  expect_equal(pack(c(13, 9, 3), 4), 
               as.raw(c(217, 48)))
  expect_equal(pack(c(11, 11, 13), 4), 
               as.raw(c(187, 208)))
  expect_equal(pack(c(12, 4, 3), 4), 
               as.raw(c(196, 48)))
  expect_equal(pack(c(12, 15, 1), 4), 
               as.raw(c(207, 16)))
  expect_equal(pack(c(10, 3, 15), 4), 
               as.raw(c(163, 240)))
  expect_equal(pack(c(1, 7, 15, 10), 4), 
               as.raw(c(23, 250)))
  expect_equal(pack(c(7, 11, 2, 14), 4), 
               as.raw(c(123, 46)))
  expect_equal(pack(c(14, 4, 5, 4), 4), 
               as.raw(c(228, 84)))
  expect_equal(pack(c(4, 12, 5, 10), 4), 
               as.raw(c(76, 90)))
  expect_equal(pack(c(1, 11, 10, 8), 4), 
               as.raw(c(27, 168)))
  expect_equal(pack(c(5, 11, 13, 14, 12, 4, 6, 8), 4), 
               as.raw(c(91, 222, 196, 104)))
  expect_equal(pack(c(6, 5, 6, 7, 2, 4, 6, 8), 4), 
               as.raw(c(101, 103, 36, 104)))
  expect_equal(pack(c(14, 2, 3, 15, 3, 13, 9, 10), 4), 
               as.raw(c(226, 63, 61, 154)))
  expect_equal(pack(c(6, 9, 4, 8, 7, 14, 2, 14, 13, 7, 5, 1, 3, 12, 8, 9, 3, 14, 11, 8, 15, 1, 1, 3, 9, 15, 5, 8, 13), 4), 
               as.raw(c(105, 72, 126, 46, 215, 81, 60, 137, 62, 184, 241, 19, 159, 88, 208)))
  expect_equal(pack(c(14, 2, 8, 13, 10, 6, 5, 10, 9, 14, 11, 3, 5), 4), 
               as.raw(c(226, 141, 166, 90, 158, 179, 80)))
  expect_equal(pack(c(3, 10, 9, 9, 3, 7, 14, 9, 6, 3, 9, 5, 11, 10, 3, 11, 6, 7, 9, 4, 4, 13, 2, 14, 9, 12), 4), 
               as.raw(c(58, 153, 55, 233, 99, 149, 186, 59, 103, 148, 77, 46, 156)))
  expect_equal(pack(c(10, 2, 11, 4, 15, 13, 6, 11, 1, 6, 2, 6, 10, 9, 1, 3, 12, 9, 2, 9, 9, 11, 6, 5), 4), 
               as.raw(c(162, 180, 253, 107, 22, 38, 169, 19, 201, 41, 155, 101)))
  expect_equal(pack(c(4, 11, 11, 8, 14, 2, 15, 15, 11, 6, 14, 11, 11, 14, 11, 10, 9, 15, 12, 3, 15, 6, 2), 4), 
               as.raw(c(75, 184, 226, 255, 182, 235, 190, 186, 159, 195, 246, 32)))
  expect_equal(pack(c(11, 4, 11, 2, 13, 7, 8, 15, 12, 9, 4), 4), 
               as.raw(c(180, 178, 215, 143, 201, 64)))
  expect_equal(pack(c(14, 8, 10, 11, 2, 12, 4, 13, 2, 13, 7, 7, 5, 5, 11, 14, 5, 2, 7, 8, 6, 15, 7, 6, 15, 9, 9, 9, 7, 15, 4, 1, 15, 3, 4, 11), 4), 
               as.raw(c(232, 171, 44, 77, 45, 119, 85, 190, 82, 120, 111, 118, 249, 153, 127, 65, 243, 75)))
  expect_equal(pack(c(7, 10, 7, 12, 10, 10, 15, 6, 12, 13, 12, 1, 12, 4), 4), 
               as.raw(c(122, 124, 170, 246, 205, 193, 196)))
  expect_equal(pack(c(9, 13, 14, 11, 7, 7, 2, 12, 4, 6, 8, 11, 12, 13, 8, 13, 8), 4), 
               as.raw(c(157, 235, 119, 44, 70, 139, 205, 141, 128)))
  expect_equal(pack(c(1, 13, 7, 5, 12, 11, 2, 7, 9, 15, 7, 14, 6, 6, 7, 9, 12, 6, 3, 8, 10, 9, 10, 3, 3, 2, 6, 13, 6, 3, 15, 11, 10, 13), 4), 
               as.raw(c(29, 117, 203, 39, 159, 126, 102, 121, 198, 56, 169, 163, 50, 109, 99, 251, 173)))
}) 

test_that("integers are unpacked correctly with alph_size = 2", {
  expect_equal(unpack(as.raw(64), 2), as.raw(c(1, 0, 0, 0)))
  expect_equal(unpack(as.raw(80), 2), as.raw(c(1, 1, 0, 0)))
  expect_equal(unpack(as.raw(84), 2), as.raw(c(1, 1, 1, 0)))
  expect_equal(unpack(as.raw(85), 2), as.raw(c(1, 1, 1, 1)))
  expect_equal(unpack(as.raw(192), 2), as.raw(c(3, 0, 0, 0)))
  expect_equal(unpack(as.raw(108), 2), as.raw(c(1, 2, 3, 0)))
  expect_equal(unpack(as.raw(255), 2), as.raw(c(3, 3, 3, 3)))
  
  expect_equal(unpack(as.raw(c(86, 104)), 2), as.raw(c(1, 1, 1, 2, 1, 2, 2, 0)))
  expect_equal(unpack(as.raw(c(122, 192)), 2), as.raw(c(1, 3, 2, 2, 3, 0, 0, 0)))
  expect_equal(unpack(as.raw(c(118, 176)), 2), as.raw(c(1, 3, 1, 2, 2, 3, 0, 0)))
  expect_equal(unpack(as.raw(c(127, 164)), 2), as.raw(c(1, 3, 3, 3, 2, 2, 1, 0)))
  expect_equal(unpack(as.raw(c(85, 156)), 2), as.raw(c(1, 1, 1, 1, 2, 1, 3, 0)))
  expect_equal(unpack(as.raw(c(110, 112)), 2), as.raw(c(1, 2, 3, 2, 1, 3, 0, 0)))
  expect_equal(unpack(as.raw(c(246, 240)), 2), as.raw(c(3, 3, 1, 2, 3, 3, 0, 0)))
  expect_equal(unpack(as.raw(c(186, 221)), 2), as.raw(c(2, 3, 2, 2, 3, 1, 3, 1)))
  expect_equal(unpack(as.raw(c(126, 124)), 2), as.raw(c(1, 3, 3, 2, 1, 3, 3, 0)))
  expect_equal(unpack(as.raw(c(187, 112)), 2), as.raw(c(2, 3, 2, 3, 1, 3, 0, 0)))
  
  expect_equal(unpack(as.raw(c(249, 238, 148)), 2), as.raw(c(3, 3, 2, 1, 3, 2, 3, 2, 2, 1, 1, 0)))
  expect_equal(unpack(as.raw(c(235, 175, 254)), 2), as.raw(c(3, 2, 2, 3, 2, 2, 3, 3, 3, 3, 3, 2)))
  expect_equal(unpack(as.raw(c(158, 187, 235)), 2), as.raw(c(2, 1, 3, 2, 2, 3, 2, 3, 3, 2, 2, 3)))
  expect_equal(unpack(as.raw(c(233, 245, 100)), 2), as.raw(c(3, 2, 2, 1, 3, 3, 1, 1, 1, 2, 1, 0)))
  expect_equal(unpack(as.raw(c(95, 221, 221)), 2), as.raw(c(1, 1, 3, 3, 3, 1, 3, 1, 3, 1, 3, 1)))
  
  expect_equal(unpack(as.raw(c(175, 223, 235, 255, 240)), 2), 
               as.raw(c(2, 2, 3, 3, 3, 1, 3, 3, 3, 2, 2, 3, 3, 3, 3, 3, 3, 3, 0, 0)))
  expect_equal(unpack(as.raw(c(247, 106, 169, 250, 186)), 2), 
               as.raw(c(3, 3, 1, 3, 1, 2, 2, 2, 2, 2, 2, 1, 3, 3, 2, 2, 2, 3, 2, 2)))
  expect_equal(unpack(as.raw(c(87, 190, 230, 170, 104)), 2), 
               as.raw(c(1, 1, 1, 3, 2, 3, 3, 2, 3, 2, 1, 2, 2, 2, 2, 2, 1, 2, 2, 0)))
})

test_that("integers are unpacked correctly with alph_size = 3", {
  expect_equal(unpack(as.raw(c(41, 128)), 3), as.raw(c(1, 2, 3, 0, 0, 0)))
  expect_equal(unpack(as.raw(c(61, 166)), 3), as.raw(c(1, 7, 3, 2, 3, 0)))
  expect_equal(unpack(as.raw(c(208, 128)), 3), as.raw(c(6, 4, 1, 0, 0, 0)))
})