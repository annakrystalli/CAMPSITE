test_that("lineage classes work", {
  
  test_lineage <- new_lineage(c(1, 0, 0, -1, 1, 0, 0))
  expect_type(test_lineage, "double")
  expect_s3_class(test_lineage, "cs_lineage")
  expect_named(test_lineage, expected = c("parent_node", "descendant_node", "start_time", "end_time", 
                                          "status", "spec_or_ext_ct", "spec_ct"))
  expect_length(test_lineage, 7)
  
  test_lineages <- new_lineages(
    test_lineage, 
    new_lineage(c(1, 0, 0, -1, -1, 0, NA))
  )
  expect_equal(nrow(test_lineages), 2)
  expect_equal(ncol(test_lineages), 7)
  expect_s3_class(test_lineages, "cs_lineages")
  expect_equal(rownames(test_lineages), c("lin1", "lin2"))

})


test_that("lineage class methods work", {

  lin_methods <- new_lineages(
    new_lineage(c(1, 2, 0, 0.1, 1, 0, 0)), 
    new_lineage(c(1, 0, 0, -1, -1, 0, NA)),
    new_lineage(c(2, 0, 0.1, -1, 1, 0.1, 0.1))
    )
  
  lin_methods <- lineage_add(lin_methods,
                             new_lineage(c(2, 0, 0.1, -1, -1, NA, NA)))
  
  expect_equal(nrow(lin_methods), 4)
  expect_equal(lineage_last_incipient(lin_methods), "lin4")
  expect_equal(lineage_last_parent(lin_methods), "lin3")
  expect_equal(lineage_next(lin_methods), 5)
  expect_equal(lineage_last_speciated(lin_methods), c("lin3", "lin4"))
  expect_equal(lineage_last_incipient(lin_methods), "lin4")
  expect_equal(lineage_last_parent(lin_methods), "lin3")
  expect_equal(lineages_active(lin_methods), c("lin2", "lin3", "lin4"))
  
  
  lin_methods[3, c(4, 5, 6)] <- c(0.2, -2, 0.2)
  lin_methods[4, c(5, 6, 7)] <- c(1, 0.2, 0.2)
  expect_equal(lineages_active(lin_methods), c("lin2", "lin4"))
  expect_equal(lineages_born(lin_methods), c("lin3", "lin4"))
  expect_equal(lineages_dead(lin_methods), c("lin1", "lin3"))
  expect_true(is_good(lin_methods, 1))
  expect_false(is_good(lin_methods, 2))
  expect_equal(lineages_edges(lin_methods), c("lin2", "lin3", "lin4"))
    
})

