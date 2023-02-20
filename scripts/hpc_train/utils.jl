function _copy_full_system_to_train_directory(
    SCRATCH_PATH,
    project_folder,
    train_folder,
    system_name,
)
    mkpath(
        joinpath(
            SCRATCH_PATH,
            project_folder,
            train_folder,
            PowerSimulationNODE.INPUT_SYSTEM_FOLDER_NAME,
        ),
    )
    cp(
        joinpath(SCRATCH_PATH, project_folder, "systems", string(system_name, ".json")),
        joinpath(
            SCRATCH_PATH,
            project_folder,
            train_folder,
            PowerSimulationNODE.INPUT_SYSTEM_FOLDER_NAME,
            string(system_name, ".json"),
        ),
        force = true,
    )
    cp(
        joinpath(
            SCRATCH_PATH,
            project_folder,
            "systems",
            string(system_name, "_validation_descriptors.json"),
        ),
        joinpath(
            SCRATCH_PATH,
            project_folder,
            train_folder,
            PowerSimulationNODE.INPUT_SYSTEM_FOLDER_NAME,
            string(system_name, "_validation_descriptors.json"),
        ),
        force = true,
    )
end
