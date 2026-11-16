*! Shared concrete frontend return type; no fit or data manipulation.
version 19.5
mata:
struct cdid_public_fit_result {
    struct cdid_result_state scalar payload
    struct cdid_bridge_request scalar inference
}
end
