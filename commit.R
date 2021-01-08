library(ggplot2)
library(data.table)
library(cowplot)
library(stringr)

# previous_data = existing
# new_data = rlang::duplicate(existing)
# new_data[between(date, "2020-07-01", "2020-07-06"), value := value + 5]
# new_data[between(date, "2020-07-01", "2020-07-06")]
# new_data = rbind(new_data, new_data[date == "2020-07-01", .(date = ymd("2020-07-15"), location, indicator, value, description)])
# new_data = new_data[date > "2020-03-15"]
#
# previous_data = rbind(previous_data, data.table(date = ymd("2020-05-05"), location = "Wales", indicator = "some_old_thing", value = 42, description = "Who cares?"))
# new_data = rbind(new_data, data.table(date = ymd("2020-05-05"), location = "East of England", indicator = "some_new_thing", value = 42, description = "Who knows?"))

# TODO
# 4. add age groups as well!?
commit = function(previous_data, new_data)
{
    # Make previous and new data into data.tables
    setDT(previous_data)
    setDT(new_data)

    # TODO check all columns expected and check no value is NA

    # Enumerate all indicators in previous and new data
    indicators = union(previous_data[, unique(indicator)], new_data[, unique(indicator)]);

    # Create final data interactively, one indicator at a time
    final_data = NULL;
    for (ind in indicators) {
        # Extract previous data for this indicator
        previous_indicator = previous_data[indicator == ind];
        setnames(previous_indicator, "value", "previous_value");

        # Extract new data for this indicator
        new_indicator = new_data[indicator == ind];
        setnames(new_indicator, "value", "new_value");

        # Merge data and mark type of change
        merged_indicator = merge(previous_indicator, new_indicator, by = c("date", "location", "indicator", "description"), all = T);
        merged_indicator[is.na(previous_value) & !is.na(new_value), change := "added"];
        merged_indicator[!is.na(previous_value) & is.na(new_value), change := "removed"];
        merged_indicator[!is.na(previous_value) & !is.na(new_value) & previous_value == new_value, change := "unchanged"];
        merged_indicator[!is.na(previous_value) & !is.na(new_value) & previous_value != new_value, change := "updated"];

        # Build plot for showing changes
        plot = ggplot();
        manual_colours = character(0);
        if (merged_indicator[change == "added", .N] > 0) {
            plot = plot +
                geom_point(data = merged_indicator[change == "added"],   aes(x = date, y = new_value,      colour = "added", shape = "added"));
            manual_colours = c(manual_colours, "#009900");
        }
        if (merged_indicator[change == "removed", .N] > 0) {
            plot = plot +
                geom_point(data = merged_indicator[change == "removed"], aes(x = date, y = previous_value, colour = "removed", shape = "removed"));
            manual_colours = c(manual_colours, "#cc0000");
        }
        if (merged_indicator[change == "unchanged", .N] > 0) {
            plot = plot +
                geom_line (data = merged_indicator[change == "unchanged"], aes(x = date, y = new_value, colour = "unchanged")) +
                geom_point(data = merged_indicator[change == "unchanged"], aes(x = date, y = new_value, colour = "unchanged", shape = "unchanged"), size = 0);
            manual_colours = c(manual_colours, "#bbbbbb");
        }
        if (merged_indicator[change == "updated", .N] > 0) {
            plot = plot +
                geom_point(data = merged_indicator[change == "updated"], aes(x = date, y = previous_value, colour = "updated (old)", shape = "updated (old)"), size = 1) +
                geom_point(data = merged_indicator[change == "updated"], aes(x = date, y = new_value,      colour = "updated (new)", shape = "updated (new)"), size = 0.5);
            manual_colours = c(manual_colours, "#6666ff", "#6666ff");
        }
        plot = plot +
            facet_wrap(~location, scales = "free_x") +
            scale_x_date(limits = merged_indicator[, c(min(date), max(date))], date_labels = "%b", date_breaks = "1 month") +
            labs(x = "date", y = ind, title = ind, subtitle = str_wrap(merged_indicator[, description[1]])) +
            theme_cowplot(font_size = 8) +
            theme(legend.position = c(0.7, 0.15), strip.background = element_blank()) +
            scale_colour_manual(values = c(added = "#009900", removed = "#cc0000", unchanged = "#bbbbbb", "updated (old)" = "#6666ff", "updated (new)" = "#6666ff"), guide = "none") +
            scale_shape_manual(values = c(added = 3, removed = 4, unchanged = 15, "updated (old)" = 1, "updated (new)" = 16),
                      guide = guide_legend(title = NULL, override.aes = list(colour = manual_colours, size = 1)))

        # Show plot and prompt user
        print(plot);
        cat(paste0("\n", "Indicator: ", ind, "\n", merged_indicator[, description[1]], "\n"));
        print(merged_indicator[, table(changes = change)]);
        a = readline(prompt = "Commit changes for this indicator? (y/n/u/a) ");

        if (a %in% c("Y", "y")) {
            setnames(new_indicator, "new_value", "value");
            final_data = rbind(final_data, new_indicator);
        } else if (a %in% c("U", "u")) {
            merged_indicator[, value := ifelse(is.na(new_value), previous_value, new_value)];
            final_data = rbind(final_data, merged_indicator[, .(date, location, indicator, value, description)]);
        } else if (a %in% c("A", "a")) {
            merged_indicator[, value := ifelse(is.na(previous_value), new_value, previous_value)];
            final_data = rbind(final_data, merged_indicator[, .(date, location, indicator, value, description)]);
        } else {
            setnames(previous_indicator, "previous_value", "value");
            final_data = rbind(final_data, previous_indicator);
        }
    }

    return (final_data)
}

most_recent = function(path, tag, format_converter)
{
    files = list.files(path, tag);
    files = files[!files %like% "^~\\$"]; # Exclude "owner lock" files created by MS Office
    dates = format_converter(str_extract(str_replace_all(files, "[ -]", ""), "[0-9]{6,12}"));
    return (file.path(path, files[which.max(dates)]))
}

read_normal = function(path, sheet)
{
    if (str_sub(path, -4) %like% "xls") {
        return (data.table(read_excel(path, sheet)));
    } else {
        return (fread(path))
    }
}

read_2_row_header = function(path, sheet)
{
    # Read header, abbreviating names
    header = data.table(read_excel(path, sheet, n_max = 1));
    empty = which(names(header) %like% "\\.\\.\\.[0-9]+");
    names(header)[empty] = names(header)[empty - 1];
    append = unname(unlist(header[1]));
    names(header) = paste0(names(header), ifelse(is.na(append), "", paste0(" ", append)));
    names(header) = make.unique(abbreviate(str_replace_all(names(header), "[^.0-9A-Za-z]", "")));

    # Read data
    data = data.table(read_excel(path, sheet, skip = 2, col_names = FALSE));
    names(data) = names(header);

    return (data)
}

read_horiz_date = function(path, sheet, range, variable.name, value.name)
{
    # Read and tame data
    data = t(data.table(read_excel(path, sheet, range, col_names = FALSE)));
    data = matrix(str_trim(c(data)), nrow = nrow(data), ncol = ncol(data));

    # Cut out summary/note rows at right side of data
    rows_keep = tail(data[, 1], -1) == as.numeric(head(data[, 1], -1)) + 1;
    rows_keep = sum(rows_keep, na.rm = TRUE) + 1;
    data = data[1:(rows_keep + 1), ];

    # Reread data as data.table and fix date column
    data[1, 1] = "date";
    data_string = paste(apply(data, 1, function(v) paste(v, collapse = "\t")), collapse = "\n");
    dt = fread(data_string);
    dt[, date := ymd("1900-01-01") + date - 2];

    dt = melt(dt, id.vars = "date");
    dt[, value := as.numeric(value)];
    setnames(dt, c("variable", "value"), c(variable.name, value.name));

    return (dt)
}

read_wales_bulk_export = function(path)
{
    files = list.files(path, "BulkExport|bulk_export");
    exports = list();
    for (f in seq_along(files)) {
        exports[[f]] = fread(file.path(path, files[f]));
    }
    exports = rbindlist(exports);
    exports = exports[!duplicated(exports, by = 1:7)];
    exports = exports[order(UpdateDateTime, HealthBoard, Hospital, Dataset, Section, Question, Measure)];
    exports[, date := as.Date(ymd_hms(UpdateDateTime))];
    return (exports)
}

melt_annotate = function(data)
{
    melted = melt(data, id.vars = c("date", "location"), variable.name = "indicator", value.name = "value");
    melted[indicator == "death_inc_line", description := "All deaths (by date of death)"];
    melted[indicator == "hospital_inc", description := "New and newly confirmed patients in hospital"];
    melted[indicator == "hospital_prev", description := "Total beds occupied"];
    melted[indicator == "icu_prev", description := "ICU beds occupied"];

    return (melted)
}

blank_fill = function(data, date_col, fill_value)
{
    date_min = data[, min(get(date_col), na.rm = T)];
    date_max = data[, max(get(date_col), na.rm = T)];
    new_data = data.table(date000 = seq(date_min, date_max, by = "1 day"));
    setnames(new_data, "date000", date_col);
    new_data = merge(new_data, data, by = date_col, all = T)
    new_data[is.na(new_data)] = fill_value;
    return (new_data)
}

