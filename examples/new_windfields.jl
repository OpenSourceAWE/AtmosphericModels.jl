using AtmosphericModels, KiteUtils

set_data_path("data") 
set = load_settings("system.yaml")
am::AtmosphericModel = AtmosphericModel(set; nowindfield=true)

new_windfields(am)
nothing
